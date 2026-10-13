import zipfile, re, sys
p = r"D:/Cander/School/3-2/system_bio/Benedicto2025_midterm_template.pptx"
z = zipfile.ZipFile(p)
names = sorted([n for n in z.namelist() if re.match(r'ppt/slides/slide\d+\.xml$', n)], key=lambda s: int(re.findall(r'\d+', s)[0]))
for n in names:
    x = z.read(n).decode('utf8')
    paras = re.findall(r'<a:p>(.*?)</a:p>', x, flags=re.S)
    out = []
    for pa in paras:
        t = ''.join(re.findall(r'<a:t>(.*?)</a:t>', pa, flags=re.S))
        if t.strip():
            out.append(t)
    print('=== ', n)
    print('\n'.join(out))
# notes
for n in sorted([n for n in z.namelist() if re.match(r'ppt/notesSlides/notesSlide\d+\.xml$', n)]):
    x = z.read(n).decode('utf8')
    t = ' '.join(re.findall(r'<a:t>(.*?)</a:t>', x, flags=re.S))
    print('--- notes', n, t[:600])
