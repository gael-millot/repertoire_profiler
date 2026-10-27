#!/usr/bin/env python3
"""
Fill a copy of the repertoire_profiler nextflow.config with user values.

The parameters of the "Parameters that must be set by the user" section of
nextflow.config are written as `name = "value" // comment` lines, inside or
outside env {} blocks. They cannot be overridden by --name on the nextflow
command line, so each web or Galaxy job gets its own filled copy of the file,
run with `nextflow -C job.config run main.nf`.

Two kinds of values:
    --values params.json   values typed by a (public, untrusted) user. Strictly
                           checked: no quote, backslash, dollar or newline, so
                           that nothing can be interpolated or executed by
                           Groovy when nextflow reads the config.
    --set name=value       values set by the server or the Galaxy wrapper
                           (trusted, e.g. "$projectDir/results").

Usage:
    fill_config.py --template nextflow.config --spec > spec.json
    fill_config.py --template nextflow.config --values params.json \
        --set system_exec=local --set out_path_ini=/tmp/res > job.config
"""

import argparse
import json
import re
import sys

PARAM_RE = re.compile(
    r"^(?P<indent>\s*)(?P<name>[A-Za-z_][A-Za-z0-9_]*)(?P<eq>\s*=\s*)"
    r"(?P<value>\"[^\"]*\"|'[^']*'|[0-9]+)(?P<rest>.*)$"
)
SECTION_RE = re.compile(r"^##\s+(?P<title>\S.*?)\s*##$")
END_USER_SECTION = "End Parameters that must be set by the user"
FORBIDDEN_CHARS = set("\"'\\$\n\r`")


def parse(lines):
    """Return the list of user parameters of nextflow.config, in file order.

    Each parameter is a dict with: name, default, quote, section, help, line.
    """
    params = []
    section = ""
    in_comment_block = False
    for i, line in enumerate(lines):
        stripped = line.strip()
        if END_USER_SECTION in stripped:
            break
        if stripped.startswith("/*"):
            in_comment_block = True
        if in_comment_block:
            m = SECTION_RE.match(stripped)
            if m and "Parameters that must be set" not in m.group("title"):
                section = m.group("title")
            if stripped.endswith("*/"):
                in_comment_block = False
            continue
        m = PARAM_RE.match(line.rstrip("\r\n"))
        if m:
            value = m.group("value")
            quote = value[0] if value[0] in "\"'" else ""
            rest = m.group("rest").strip()
            help_text = rest[2:].strip() if rest.startswith("//") else ""
            params.append({
                "name": m.group("name"),
                "default": value[1:-1] if quote else value,
                "quote": quote,
                "section": section,
                "help": help_text,
                "line": i,
            })
        elif stripped.startswith("//") and params and params[-1]["line"] < i and not stripped.startswith("///"):
            # continuation lines of the previous parameter comment (e.g. align_seq), only inside env {}
            prev = lines[params[-1]["line"] + 1:i]
            if all(p.strip().startswith("//") for p in prev) and line.startswith((" ", "\t")):
                params[-1]["help"] += "\n" + stripped[2:].strip()
    return params


def check_user_value(name, value):
    if not isinstance(value, (str, int, float)) or isinstance(value, bool):
        raise ValueError(f"{name}: value must be a string or a number")
    value = str(value)
    bad = FORBIDDEN_CHARS.intersection(value)
    if bad:
        raise ValueError(f"{name}: forbidden character(s) {''.join(sorted(bad))!r}")
    if len(value) > 1000:
        raise ValueError(f"{name}: value too long")
    return value


def fill(lines, params, values):
    """Return the lines with the values replaced. values: {name: str}"""
    by_name = {p["name"]: p for p in params}
    unknown = set(values) - set(by_name)
    if unknown:
        raise ValueError("unknown parameter(s): " + ", ".join(sorted(unknown)))
    out = list(lines)
    for name, value in values.items():
        p = by_name[name]
        if p["quote"]:
            if p["quote"] in value or "\n" in value:
                raise ValueError(f"{name}: value cannot contain {p['quote']} or a newline")
            new_value = p["quote"] + value + p["quote"]
        elif re.fullmatch(r"[0-9]+", value):
            new_value = value
        else:
            raise ValueError(f"{name}: integer expected")
        line = out[p["line"]]
        eol = line[len(line.rstrip("\r\n")):]
        m = PARAM_RE.match(line.rstrip("\r\n"))
        out[p["line"]] = m.group("indent") + name + m.group("eq") + new_value + m.group("rest") + eol
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--template", required=True, help="path of the nextflow.config of repertoire_profiler")
    parser.add_argument("--spec", action="store_true", help="print the user parameters as JSON and exit")
    parser.add_argument("--values", help="JSON file of untrusted user values {name: value}")
    parser.add_argument("--set", action="append", default=[], metavar="NAME=VALUE", help="trusted value, can be repeated")
    parser.add_argument("--locked", help="file of parameter names that users cannot set (one per line, # for comments)")
    parser.add_argument("-o", "--output", help="output file (default: stdout)")
    args = parser.parse_args()

    with open(args.template, encoding="utf-8", newline="") as f:
        lines = f.readlines()
    params = parse(lines)

    if args.spec:
        json.dump(params, sys.stdout, indent=1)
        sys.stdout.write("\n")
        return

    locked = set()
    if args.locked:
        with open(args.locked, encoding="utf-8") as f:
            locked = {l.split("#")[0].strip() for l in f} - {""}

    values = {}
    try:
        if args.values:
            with open(args.values, encoding="utf-8") as f:
                user_values = json.load(f)
            if not isinstance(user_values, dict):
                raise ValueError("the --values file must contain a JSON object")
            forbidden = locked.intersection(user_values)
            if forbidden:
                raise ValueError("locked parameter(s): " + ", ".join(sorted(forbidden)))
            values = {k: check_user_value(k, v) for k, v in user_values.items()}
        for item in args.set:
            name, sep, value = item.partition("=")
            if not sep:
                raise ValueError(f"--set {item}: NAME=VALUE expected")
            values[name] = value
        out = fill(lines, params, values)
    except ValueError as e:
        sys.exit(f"fill_config.py: ERROR: {e}")

    if args.output:
        with open(args.output, "w", encoding="utf-8", newline="") as f:
            f.writelines(out)
    else:
        sys.stdout.writelines(out)


if __name__ == "__main__":
    main()
