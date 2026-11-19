# Contributing to ApproxKit

Thanks for your interest. ApproxKit provides bit-accurate Python models, error metrics, a Pareto explorer and Verilog for approximate adders and multipliers.

## Getting started

```
pip install -e ".[plots,test]"        # numpy; matplotlib for plots; pytest optional
python examples/run_dse.py --n 8
python examples/demo_apps.py
```

## Testing

```
python tests/test_units.py            # plain asserts (pytest tests also works)
python scripts/rtl_regress.py --out build --synth   # Python models vs Verilog; needs iverilog and yosys
```

The Python model and the Verilog must stay bit-exact; the RTL regression is the check for that.

## Pull requests

- Open an issue first for larger changes.
- Keep PRs focused and fill in the pull request template.
- Follow `.editorconfig` and match the style of surrounding code.
- Say which simulator and tool versions you used for RTL changes.

## Conduct

Participation is governed by the [Code of Conduct](CODE_OF_CONDUCT.md).
