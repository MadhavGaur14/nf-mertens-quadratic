"""Static checks for main.tex (no LaTeX needed): environments, citations, references, braces, CHECK marks.
Prints the counts so that an empty match cannot pass silently."""
import collections, re, sys
from pathlib import Path

t = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).with_name('main.tex')).read_text(encoding='utf-8')
b = collections.Counter(re.findall(r'\\begin\{(\w+\*?)\}', t))
e = collections.Counter(re.findall(r'\\end\{(\w+\*?)\}', t))
print(f"environments: {sum(b.values())} begins, {sum(e.values())} ends; unbalanced:",
      {k: (b[k], e[k]) for k in set(b) | set(e) if b[k] != e[k]} or 'none')
keys = set(re.findall(r'\\bibitem\[[^\]]*\]\{([^}]+)\}', t))
cites = {c.strip() for g in re.findall(r'\\cite(?:\[[^\]]*\])?\{([^}]+)\}', t) for c in g.split(',')}
print(f"bibitems: {len(keys)}, distinct cited keys: {len(cites)}; cited but missing:", cites - keys or 'none',
      "| uncited:", keys - cites or 'none')
labels = set(re.findall(r'\\label\{([^}]+)\}', t))
refs = set(re.findall(r'\\(?:ref|eqref)\{([^}]+)\}', t))
print(f"labels: {len(labels)}, distinct refs: {len(refs)}; refs without label:", refs - labels or 'none')
body = re.sub(r'(?<!\\)%.*', '', t)
print("brace balance (should be 0):", body.count('{') - body.count('}'))
print("red CHECK marks remaining:", len(re.findall(r'\\checkme\{', t)))
