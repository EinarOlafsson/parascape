"""Bounded local raw cache with checked, temporary NAS overflow.

Retains the original cache receipt/lock API used by the running RNA worker.
Owned symlinks expose checked NAS archives at their original local names. No
quantification settings or output values change. Raw overflow is removed only
after the independent queue verifier accepts every planned use of a run.
"""
import argparse
from collections import defaultdict
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import uuid
import omics_download_cache as L
import omics_quant as Q

GIB=1024**3
OWNER='ParaScape task-owned temporary verified raw-download overflow v1\n'


def regular_bytes(root):
    return sum(p.stat().st_size for p in root.iterdir()
        if not p.is_symlink() and p.is_file() and p.suffix in ('.sra','.part'))


def _runs(active):
    """Active runs as a set: one per quantifier worker; accepts a single accession or None."""
    return set(active) if isinstance(active, (set, frozenset, list, tuple)) else ({active} if active else set())


def nas_cap_blocks(stored, incoming, cap):
    """True when adding `incoming` bytes would take the NAS tier past its cap (None = uncapped).

    The NAS is shared with other users; the free-space reserve alone let this tier grow
    from 214 GB to 1.58 TB in one night, so a fixed byte budget is enforced separately.
    """
    return cap is not None and stored + incoming > cap


def owned_pointer(path, target):
    return path.is_symlink() and os.path.abspath(os.readlink(path))==str(target.absolute())


class Cache:
    def __init__(self,local,nas,reserve,max_bytes,events):
        self.local=Path(local);self.nas=Path(nas);self.reserve=reserve;self.max_bytes=max_bytes;self.events=Path(events)

    def event(self,operation,**fields):
        row=dict(epoch=time.time(),operation=operation,**fields)
        self.events.parent.mkdir(parents=True,exist_ok=True)
        with self.events.open('a') as stream:stream.write(json.dumps(row)+'\n')

    def pair(self,run):
        if not re.fullmatch(r'[SED]RR\d+',run):raise ValueError('Unsafe run ID')
        return self.local/(run+'.sra'),self.nas/(run+'.sra')

    def receipt(self,path,url=None,size=None):
        receipt=json.loads(path.with_suffix('.json').read_bytes())
        if (not isinstance(receipt['bytes'],int) or receipt['bytes']<=0
                or not re.fullmatch('[0-9a-f]{64}',receipt['sha256'])
                or (url is not None and receipt['url']!=url)
                or (size is not None and receipt['bytes']!=size)):
            raise ValueError('Raw archive receipt differs')
        if not path.is_file() or path.stat().st_size!=receipt['bytes']:
            raise ValueError('Raw archive size differs')
        return receipt

    def publish(self,run,receipt):
        local,nas=self.pair(run)
        if nas.is_symlink():raise ValueError('NAS payload must be ordinary')
        self.receipt(nas,receipt['url'],receipt['bytes'])
        if local.is_symlink() and not owned_pointer(local,nas):raise ValueError('Unowned raw pointer')
        L.save(local.with_suffix('.json'),receipt)
        temporary=local.with_name(local.name+'.'+uuid.uuid4().hex+'.link')
        temporary.symlink_to(nas.absolute());temporary.replace(local)

    def spill(self,run,active=None):
        """Caller can never evict the active run; verify destination before unlink."""
        if run in _runs(active):return False
        local,nas=self.pair(run)
        with L.archive_lock(local):
            if local.is_symlink():
                if not owned_pointer(local,nas):raise ValueError('Unowned cache symlink')
                return False
            if not local.exists():return False
            receipt=self.receipt(local)
            signature=(local.stat().st_size,local.stat().st_mtime_ns,local.stat().st_ctime_ns)
            if nas.exists():
                # A re-download of a run already on the NAS gets a new 'finished' time; content must match.
                same=lambda r:{k:r.get(k) for k in ('url','bytes','sha256')}
                if nas.is_symlink() or same(self.receipt(nas))!=same(receipt) or L.digest(nas)!=receipt['sha256']:
                    raise ValueError('Existing NAS raw object differs')
            else:
                if shutil.disk_usage(self.nas).free-receipt['bytes']<self.reserve:
                    raise OSError('NAS reserve prevents spill')
                temporary=nas.with_name(nas.name+'.moving');digest=hashlib.sha256();copied=0
                with local.open('rb') as source,temporary.open('wb') as target:
                    for block in iter(lambda:source.read(1024**2),b''):
                        if shutil.disk_usage(self.nas).free-len(block)<self.reserve:
                            raise OSError('NAS reserve reached during spill')
                        target.write(block);digest.update(block);copied+=len(block)
                    target.flush();os.fsync(target.fileno())
                if copied!=receipt['bytes'] or digest.hexdigest()!=receipt['sha256'] or L.digest(temporary)!=receipt['sha256']:
                    raise ValueError('Copied NAS raw checksum differs; local retained')
                if signature!=(local.stat().st_size,local.stat().st_mtime_ns,local.stat().st_ctime_ns):
                    raise ValueError('Local raw archive changed during spill')
                temporary.replace(nas);L.save(nas.with_suffix('.json'),receipt)
            self.publish(run,receipt)
            self.event('spill_verified_local_raw_to_NAS',run=run,source=str(local),destination=str(nas),
                bytes=receipt['bytes'],sha256=receipt['sha256'],destination_independently_rehashed=True)
            return True

    def cleanup(self,plan,verified,active):
        uses=defaultdict(set)
        for row in plan:uses[row['run_accession']].add(row['study']+'/'+row['run_accession'])
        for run,required in uses.items():
            if run in _runs(active) or not required.issubset(verified):continue
            local,nas=self.pair(run)
            with L.archive_lock(local):
                if local.is_symlink() and not owned_pointer(local,nas):raise ValueError('Unowned raw pointer')
                if nas.exists():
                    receipt=self.receipt(nas)
                    if local.exists() and not owned_pointer(local,nas):
                        continue  # Never delete an unrelated/replaced primary object.
                    if owned_pointer(local,nas):local.unlink()
                    nas.unlink();nas.with_suffix('.json').unlink()
                    local.with_suffix('.json').unlink(missing_ok=True)
                    self.event('remove_temporary_NAS_raw_after_all_planned_uses_verified',run=run,
                        bytes=receipt['bytes'],sha256=receipt['sha256'],verified_study_run_keys=sorted(required))
                elif local.is_file() and not local.is_symlink():
                    receipt=self.receipt(local);local.unlink();local.with_suffix('.json').unlink()
                    self.event('remove_completed_local_raw',run=run,sha256=receipt['sha256'])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,required=True)
    parser.add_argument('--nas',type=Path,required=True)
    parser.add_argument('--max-gib',type=int,default=200)
    parser.add_argument('--reserve-gib',type=int,default=25)
    # The NAS is shared; without a cap the overflow tier grows until only the reserve is left.
    parser.add_argument('--nas-max-gib',type=int,default=0,help='hard cap on raw bytes kept in the NAS tier (0 = no cap)')
    # Local scratch (mnt1) is shared with tutorial recordings; with room on the NAS, download there first.
    parser.add_argument('--prefer-nas',action='store_true',help='download to the NAS tier first; local only when the NAS cap is reached')
    args=parser.parse_args()
    if min(args.max_gib,args.reserve_gib)<=0:raise ValueError('Invalid raw storage bounds')
    local=args.cache.resolve(strict=True);nas=args.nas
    filesystem=subprocess.check_output(['findmnt','-T',str(nas.parent),'-n','-o','FSTYPE'],text=True).strip()
    if filesystem not in ('nfs','nfs4','cifs'):raise OSError('NAS is not mounted as network storage')
    if nas.is_symlink():raise ValueError('Linked NAS overflow root')
    nas.mkdir(exist_ok=True)
    marker=nas/'owner.txt'
    if marker.exists():
        if marker.read_text()!=OWNER:raise ValueError('NAS cache ownership differs')
    elif any(nas.iterdir()):raise ValueError('Nonempty unowned NAS cache')
    else:marker.write_text(OWNER)
    with (local/'.prefetch.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        plan=[json.loads(line) for line in Path(Q.PLAN).read_bytes().splitlines() if line.strip()]
        for row in plan:
            if not re.fullmatch('[A-Za-z0-9_-]+',row['study']):raise ValueError('Unsafe study')
        cache=Cache(local,nas,args.reserve_gib*GIB,args.max_gib*GIB,nas/'operations.jsonl')
        # The legacy consumer follows the published pointers unchanged. Its
        # producer budget must count only bytes physically stored in each tier.
        L.cache_bytes=regular_bytes
        for row in plan:cache.pair(row['run_accession'])
        cache.event('NAS_overflow_worker_started',pid=os.getpid(),source_sha256=L.digest(Path(__file__)),
            local_cache=str(local),nas_cache=str(nas),local_limit_bytes=cache.max_bytes,reserve_bytes=cache.reserve,
            NAS_budget=(f'Hard cap {args.nas_max_gib} GiB of raw bytes in the NAS tier' if args.nas_max_gib else 'Available NAS filesystem space above reserve; no fixed extra cache allocation'))
        nas_cap=args.nas_max_gib*GIB if args.nas_max_gib else None
        failures={}

        def report(**fields):
            value=dict(pid=os.getpid(),max_bytes=cache.max_bytes,reserve_bytes=cache.reserve,
                rate_bytes_per_second=0,updated=time.strftime('%Y-%m-%dT%H:%M:%S%z'),
                local_cache_bytes=regular_bytes(local),nas_cache_bytes=regular_bytes(nas),**fields)
            L.save(local/'status.json',value)

        while True:
            try:
                # One active file per quantifier worker (active.json or active-w<k>.json).
                active=set()
                for path in local.glob('active*.json'):
                    try:active.add(json.loads(path.read_bytes())['run'])
                    except (OSError,KeyError,ValueError):continue
            except OSError:active=set()
            verified=json.loads((local/'automation/verified.json').read_bytes())
            cache.cleanup(plan,verified,active)
            pending=[r for r in plan if r['study']+'/'+r['run_accession'] not in verified]
            if not pending:
                report(phase='complete',reason='Every selected planned use independently verified');return
            # Leave enough local scratch room for the next cached run's existing
            # eight-times-archive expansion estimate, plus the established reserve.
            next_size=next((cache.pair(r['run_accession'])[0].stat().st_size for r in pending
                if r['run_accession'] not in active and cache.pair(r['run_accession'])[0].is_file()),0)
            scratch_required=max(cache.reserve,8*next_size+25*10**9)
            low_scratch=shutil.disk_usage(local).free<scratch_required
            if low_scratch:
                candidates=sorted((p for p in local.glob('[SED]RR*.sra') if not p.is_symlink() and p.stem not in active),
                    key=lambda p:p.stat().st_size,reverse=True)
                if candidates and nas_cap_blocks(regular_bytes(nas),candidates[0].stat().st_size,nas_cap):
                    report(phase='waiting',reason='NAS cache cap reached; local scratch is low until verified runs are cleaned up');time.sleep(300);continue
                if candidates:
                    report(phase='spilling_checked_raw_to_NAS',reason='Preserve next-run scratch expansion space',run=candidates[0].stem)
                    cache.spill(candidates[0].stem,active)
                    continue
            attempted=False
            for row in pending:
                run=row['run_accession'];archive,nas_archive=cache.pair(run)
                if run in active or failures.get(run,0)>time.monotonic():continue
                if archive.exists():
                    if archive.is_symlink() and not owned_pointer(archive,nas_archive):raise ValueError('Unowned raw cache link')
                    continue
                reason=L.pressure_reason()
                if reason:
                    report(phase='waiting',reason=reason);break
                url=Q.SRA_ODP.format(run=run);size=Q._head_size(url)
                if not size:
                    failures[run]=time.monotonic()+3600;continue
                local_ok=(regular_bytes(local)+size<=cache.max_bytes
                    and shutil.disk_usage(local).free-size>=scratch_required)
                if args.prefer_nas and not nas_cap_blocks(regular_bytes(nas),size,nas_cap):
                    local_ok=False
                target=archive if local_ok else nas_archive
                reserve=cache.reserve if target==nas_archive else scratch_required
                part=target.with_suffix('.part')
                remaining=0 if target.exists() else size-(part.stat().st_size if part.exists() else 0)
                if target==nas_archive and nas_cap_blocks(regular_bytes(nas),remaining,nas_cap):
                    report(phase='waiting',reason='NAS cache cap reached',run=run,expected_bytes=size,downloaded_bytes=0,storage_tier='NAS');break
                if shutil.disk_usage(target.parent).free-remaining<reserve:
                    report(phase='waiting',reason='NAS reserve or local scratch reserve',run=run,expected_bytes=size,
                        downloaded_bytes=0,storage_tier='NAS' if target==nas_archive else 'local');break
                report(phase='downloading',reason='',run=run,expected_bytes=size,downloaded_bytes=0,
                    storage_tier='NAS' if target==nas_archive else 'local')
                try:
                    with L.archive_lock(archive):
                        if archive.exists():continue
                        budget=cache.max_bytes if target==archive else min(regular_bytes(nas)+shutil.disk_usage(nas).free-cache.reserve,nas_cap or float('inf'))
                        if target==nas_archive and target.exists():
                            receipt=cache.receipt(target,url,size)
                            if L.digest(target)!=receipt['sha256']:raise ValueError('Recovered NAS raw checksum differs')
                        elif target==archive:
                            # Already holding the primary lock; avoid taking it twice.
                            receipt=L._download(url,target,size,rate=0,reserve=reserve,budget=budget,
                                guard=L.pressure_reason,report=lambda **kw:report(run=run,storage_tier='local',expected_bytes=size,**kw))
                        else:
                            receipt=L.download(url,target,size,rate=0,reserve=reserve,budget=budget,
                                report=lambda **kw:report(run=run,storage_tier='NAS',expected_bytes=size,**kw))
                        if target==nas_archive:cache.publish(run,receipt)
                        cache.event('download_raw_archive_complete',run=run,path=str(target),bytes=size,
                            sha256=receipt['sha256'],storage_tier='NAS' if target==nas_archive else 'local')
                    report(phase='ready',reason='',run=run,expected_bytes=size,downloaded_bytes=size,
                        storage_tier='NAS' if target==nas_archive else 'local')
                except Exception as error:
                    failures[run]=time.monotonic()+300
                    report(phase='retry_later',reason=str(error),run=run,expected_bytes=size)
                    cache.event('retry_raw_download',run=run,error=str(error))
                attempted=True;break
            if not attempted:time.sleep(30)


if __name__=='__main__':main()
