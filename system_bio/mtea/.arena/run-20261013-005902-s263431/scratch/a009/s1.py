import openpyxl, sys
sys.stdout.reconfigure(encoding='utf-8')
for f in ['MOESM2','MOESM3','MOESM4','MOESM5','MOESM7']:
    wb=openpyxl.load_workbook(f'si/{f}.xlsx',read_only=True)
    print('==',f)
    for ws in wb.worksheets:
        rows=list(ws.iter_rows(values_only=True))
        print(' sheet',ws.title,len(rows),'rows; first rows:')
        for r in rows[:4]: print('   ',r[:10])
