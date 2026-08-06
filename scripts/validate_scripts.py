"""Fast, license-free structural checks for the simulation workflow."""

from __future__ import annotations

from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
REPO = SCRIPT_DIR.parent

EXPECTED_CALLS = {
    "final_device.lsf": ("build_planar_double_junction", "prepare_generation_group"),
    "final_embedded_5nm_effective.lsf": (
        "build_embedded_active",
        "prepare_generation_group",
        "QD5nm_20pct_in_MAPbI3_MG",
    ),
    "qd_patterns/effective_medium_blend.lsf": ("build_embedded_active",),
    "qd_patterns/explicit_5nm_patterns.lsf": (
        "build_embedded_active",
        "addsphere",
        '"square"',
        '"hex"',
        '"honeycomb"',
    ),
}

REQUIRED_FILES = (
    "models/fdtd/final_planar_double_junction.fsp",
    "models/fdtd/final_embedded_5nm_effective.fsp",
    "models/charge/best_0p80eV.ldev",
    "models/charge/hedge_1p03eV.ldev",
    "data/generation/planar_fdtd_raw.mat",
    "data/generation/planar_hedge_fdtd_raw.mat",
    "data/generation/charge_best.mat",
    "data/generation/charge_hedge.mat",
    "data/generation/embedded_5nm_fdtd_raw.mat",
    "data/generation/charge_embedded.mat",
)

FORBIDDEN_TEXT = {
    'getresult("Pabs_thermal"': "thermal power cannot be carrier generation",
    "scratch_runs": "persistent scratch directories are not allowed",
    "sinusoidal": "removed sinusoidal workflow is still referenced",
    "champion": "use 'best' naming",
}


def check_lsf_balance(path: Path, text: str) -> list[str]:
    """Check delimiters while ignoring quoted strings and # comments."""
    errors: list[str] = []
    pairs = {")": "(", "]": "[", "}": "{"}
    stack: list[tuple[str, int]] = []
    quote = False
    line = 1
    index = 0
    while index < len(text):
        char = text[index]
        if char == "\n":
            line += 1
            index += 1
            continue
        if not quote and char == "#":
            newline = text.find("\n", index)
            if newline < 0:
                break
            index = newline
            continue
        if char == '"':
            quote = not quote
            index += 1
            continue
        if not quote:
            if char in "([{":
                stack.append((char, line))
            elif char in pairs:
                if not stack or stack[-1][0] != pairs[char]:
                    errors.append(f"{path.relative_to(REPO)}:{line}: unmatched {char}")
                else:
                    stack.pop()
        index += 1

    if quote:
        errors.append(f"{path.relative_to(REPO)}: unterminated quoted string")
    for char, opening_line in stack:
        errors.append(
            f"{path.relative_to(REPO)}:{opening_line}: unclosed {char}"
        )
    return errors


def main() -> None:
    errors: list[str] = []
    lsf_paths = sorted(SCRIPT_DIR.rglob("*.lsf"))
    python_paths = sorted(SCRIPT_DIR.rglob("*.py"))

    for path in (*lsf_paths, *python_paths):
        text = path.read_text(encoding="utf-8")
        if path.suffix == ".lsf":
            errors.extend(check_lsf_balance(path, text))
        if any(line.startswith("+") for line in text.splitlines()):
            errors.append(f"{path.relative_to(REPO)}: patch-artifact '+' line")
        if path.resolve() != Path(__file__).resolve():
            for token, explanation in FORBIDDEN_TEXT.items():
                if token.lower() in text.lower():
                    errors.append(f"{path.relative_to(REPO)}: {explanation}")

    for filename, calls in EXPECTED_CALLS.items():
        path = SCRIPT_DIR / filename
        if not path.is_file():
            errors.append(f"missing required script: scripts/{filename}")
            continue
        text = path.read_text(encoding="utf-8")
        for call in calls:
            if call not in text:
                errors.append(f"scripts/{filename}: missing {call}")

    for relative in REQUIRED_FILES:
        if not (REPO / relative).is_file():
            errors.append(f"missing required file: {relative}")

    for path in python_paths:
        try:
            compile(path.read_text(encoding="utf-8"), str(path), "exec")
        except SyntaxError as exc:
            errors.append(f"{path.relative_to(REPO)}:{exc.lineno}: {exc.msg}")

    if errors:
        print("Script validation: FAIL")
        for error in errors:
            print(f"  - {error}")
        raise SystemExit(1)

    print(
        f"Script validation: PASS "
        f"({len(lsf_paths)} LSF, {len(python_paths)} Python)"
    )


if __name__ == "__main__":
    main()
