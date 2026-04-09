"""
ExomeFlow — Requirements Checker

Checks all system tools and Python packages required by ExomeFlow.
Run before starting the pipeline to catch missing dependencies early.

Usage:
    python check_requirements.py
"""

from __future__ import annotations

import importlib
import shutil
import subprocess
import sys
from dataclasses import dataclass
from typing import Optional


# ── ANSI colours (fallback to plain text if not supported) ────────────────────
GREEN  = "\033[92m"
RED    = "\033[91m"
YELLOW = "\033[93m"
BOLD   = "\033[1m"
RESET  = "\033[0m"

def ok(msg: str)   -> str: return f"{GREEN}  ✔  {msg}{RESET}"
def fail(msg: str) -> str: return f"{RED}  ✘  {msg}{RESET}"
def warn(msg: str) -> str: return f"{YELLOW}  ⚠  {msg}{RESET}"
def bold(msg: str) -> str: return f"{BOLD}{msg}{RESET}"


# ── REQUIREMENT DEFINITIONS ───────────────────────────────────────────────────

@dataclass
class SystemTool:
    name: str
    min_version: str
    version_cmd: list[str]          # command to get version string
    version_pattern: str            # substring to extract from output
    install_hint: str


@dataclass
class PythonPackage:
    import_name: str
    pypi_name: str
    min_version: str


SYSTEM_TOOLS: list[SystemTool] = [
    SystemTool(
        name="bwa",
        min_version="0.7.17",
        version_cmd=["bwa"],
        version_pattern="Version:",
        install_hint="conda install -c bioconda bwa",
    ),
    SystemTool(
        name="samtools",
        min_version="1.13",
        version_cmd=["samtools", "--version"],
        version_pattern="samtools",
        install_hint="conda install -c bioconda samtools",
    ),
    SystemTool(
        name="gatk",
        min_version="4.6.0",
        version_cmd=["gatk", "--version"],
        version_pattern="",
        install_hint=(
            "Download from https://github.com/broadinstitute/gatk/releases\n"
            "    then add to PATH:  export PATH=/path/to/gatk-4.6.x.x:$PATH"
        ),
    ),
    SystemTool(
        name="fastp",
        min_version="0.20.1",
        version_cmd=["fastp", "--version"],
        version_pattern="fastp",
        install_hint="conda install -c bioconda fastp",
    ),
    SystemTool(
        name="perl",
        min_version="5.26",
        version_cmd=["perl", "--version"],
        version_pattern="v5.",
        install_hint="conda install perl   OR   sudo apt install perl",
    ),
]

PYTHON_PACKAGES: list[PythonPackage] = [
    PythonPackage("typer",      "typer",      "0.12.0"),
    PythonPackage("rich",       "rich",       "13.0.0"),
    PythonPackage("pandas",     "pandas",     "2.0.0"),
    PythonPackage("matplotlib", "matplotlib", "3.7.0"),
]


# ── VERSION HELPERS ───────────────────────────────────────────────────────────

def parse_version(text: str) -> tuple[int, ...]:
    """Extract the first dotted-integer version string found in *text*."""
    import re
    m = re.search(r"(\d+)\.(\d+)(?:\.(\d+))?", text)
    if not m:
        return (0,)
    return tuple(int(g) for g in m.groups() if g is not None)


def version_ok(found: tuple[int, ...], minimum: str) -> bool:
    min_tuple = parse_version(minimum)
    return found >= min_tuple


# ── CHECKERS ─────────────────────────────────────────────────────────────────

def check_system_tool(tool: SystemTool) -> tuple[bool, str]:
    """
    Returns (passed, message).
    """
    # 1. Is it on PATH?
    if not shutil.which(tool.name):
        return False, (
            f"{tool.name} not found on PATH\n"
            f"    Install : {tool.install_hint}"
        )

    # 2. Can we get a version?
    try:
        result = subprocess.run(
            tool.version_cmd,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=10,
        )
        output = result.stdout + result.stderr
    except Exception as exc:
        return False, f"{tool.name} found but version check failed: {exc}"

    found = parse_version(output)
    if not version_ok(found, tool.min_version):
        found_str = ".".join(str(x) for x in found)
        return False, (
            f"{tool.name} version {found_str} is below minimum {tool.min_version}\n"
            f"    Update  : {tool.install_hint}"
        )

    found_str = ".".join(str(x) for x in found)
    return True, f"{tool.name} {found_str}  (required >= {tool.min_version})"


def check_python_package(pkg: PythonPackage) -> tuple[bool, str]:
    """
    Returns (passed, message).
    """
    try:
        mod = importlib.import_module(pkg.import_name)
    except ImportError:
        return False, (
            f"{pkg.pypi_name} is not installed\n"
            f"    Install : pip install {pkg.pypi_name}>={pkg.min_version}"
        )

    # Try to get __version__
    version_str = getattr(mod, "__version__", None)
    if version_str is None:
        # Some packages use importlib.metadata
        try:
            from importlib.metadata import version
            version_str = version(pkg.pypi_name)
        except Exception:
            pass

    if version_str:
        found = parse_version(version_str)
        if not version_ok(found, pkg.min_version):
            return False, (
                f"{pkg.pypi_name} {version_str} is below minimum {pkg.min_version}\n"
                f"    Update  : pip install --upgrade {pkg.pypi_name}>={pkg.min_version}"
            )
        return True, f"{pkg.pypi_name} {version_str}  (required >= {pkg.min_version})"

    return True, f"{pkg.pypi_name} installed  (could not determine version)"


def check_annovar(hint_only: bool = True) -> tuple[bool, str]:
    """
    ANNOVAR is not on a standard PATH — check common locations and warn.
    """
    common_paths = [
        "/media/drprabudh/m1/annovar/table_annovar.pl",
        "/usr/local/bin/annovar/table_annovar.pl",
        "/opt/annovar/table_annovar.pl",
    ]
    for p in common_paths:
        import os
        if os.path.isfile(p):
            return True, f"ANNOVAR found at {p}"

    return False, (
        "ANNOVAR (table_annovar.pl) not found in common locations\n"
        "    Download: https://annovar.openbioinformatics.org\n"
        "    Then pass to CLI: --annovar-bin /path/to/annovar"
    )


def check_python_version() -> tuple[bool, str]:
    major, minor = sys.version_info[:2]
    v = f"{major}.{minor}.{sys.version_info.micro}"
    if (major, minor) < (3, 9):
        return False, (
            f"Python {v} is below minimum 3.9\n"
            f"    Update  : conda install python>=3.9"
        )
    return True, f"Python {v}  (required >= 3.9)"


# ── MAIN ─────────────────────────────────────────────────────────────────────

def main() -> None:
    passed_all = True

    print()
    print(bold("=" * 58))
    print(bold("  ExomeFlow — Requirements Check"))
    print(bold("=" * 58))

    # ── Python interpreter
    print(f"\n{bold('Python Interpreter')}")
    p, msg = check_python_version()
    print(ok(msg) if p else fail(msg))
    if not p:
        passed_all = False

    # ── Python packages
    print(f"\n{bold('Python Packages')}")
    for pkg in PYTHON_PACKAGES:
        p, msg = check_python_package(pkg)
        print(ok(msg) if p else fail(msg))
        if not p:
            passed_all = False

    # ── System tools
    print(f"\n{bold('System Tools')}")
    for tool in SYSTEM_TOOLS:
        p, msg = check_system_tool(tool)
        print(ok(msg) if p else fail(msg))
        if not p:
            passed_all = False

    # ── ANNOVAR (special case)
    p, msg = check_annovar()
    print(warn(msg) if not p else ok(msg))
    # ANNOVAR is a warning, not a hard fail (path passed via CLI)

    # ── Summary
    print()
    print(bold("=" * 58))
    if passed_all:
        print(ok(bold("All requirements satisfied — ExomeFlow is ready to run.")))
    else:
        print(fail(bold("Some requirements are missing. Fix the issues above and re-run.")))
    print(bold("=" * 58))
    print()

    sys.exit(0 if passed_all else 1)


if __name__ == "__main__":
    main()
