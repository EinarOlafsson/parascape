# 96 · The ParaScape module mould, as a skill

**Status: built 2026-11-06.** Einar agreed the Map design ("sounds good") and it is rebuilt to it
(`parascape/map_module.py`, `parascape/map_reading.py`). The mould is the skill
`~/.claude/skills/parascape-module/SKILL.md`, versioned here as `instructions/skills/parascape-module.md`
(`.claude/` is git-ignored). Next: the Gene module through the skill. Einar: "let's get this
one right then use it as a model for the rest so we can iterate faster ... explain to me what the goal of the
module is and I'll tell you what needs to change to make it more intuitive, then you can use what we come up
with as a mould for the other modules in ParaScape. So it should become a skill, parascape module skill."

## The plan
1. Review Map with Einar: its goal, and what makes it hard to understand.
2. Rebuild Map to what is agreed.
3. Write the agreed pattern down as a Claude Code skill, `parascape-module`: the goal every module serves
   (gene/list information accessible and synthesised; inference from labels or lists to knowledge), the
   shell layout, what goes left / in the preview / in actions, the wording rules, the Navigator downstream
   entry, the tests and the render check -- so each further module is built the same way, quickly.
4. Move the remaining modules with it (Gene first, synthesis first, as Einar chose).

## Gate
Einar says Map is right; the skill builds the next module to the same standard without re-explaining.
