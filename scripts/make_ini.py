"""Write a mendel-go input file from the repository defaults with key = value overrides.
Usage: python3 -I scripts/make_ini.py <defaults.ini> <out.ini> key=value [key=value ...]
Values are written verbatim (quote strings yourself, e.g. 'selection_model="spps"')."""
import sys, re
defaults, out, *overrides = sys.argv[1:]
text = open(defaults).read()
for kv in overrides:
    key, val = kv.split('=', 1)
    pat = re.compile(r'^(\s*' + re.escape(key) + r'\s*=\s*)([^#\n]*)(.*)$', re.M)
    if not pat.search(text): sys.exit(f'parameter {key!r} not found in {defaults}')
    text = pat.sub(lambda m: m.group(1) + val + '   ' + m.group(3).strip(), text, count=1)
open(out, 'w').write(text)
print(out)
