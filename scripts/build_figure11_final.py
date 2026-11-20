"""Figure 11 (storyboard 2026-11-20): host-facing effectors are the most variable and least ordered
proteins, and their host response is largely shared. Runs the effector builder (panels A-C) and
the host-response builder (D-F) in one process and composes once."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ["PP_NO_COMPOSE"] = "1"
import build_figure08_v2 as E   # noqa: E402
import build_figure10_v2 as H   # noqa: E402
import paper_panels as PP       # noqa: E402

E.main()
H.main()
os.environ.pop("PP_NO_COMPOSE")
print(PP.compose("figure_11__final", "main", "Host-facing effectors are the most variable and "
                 "least ordered proteins, and the host response to them is largely shared."))
