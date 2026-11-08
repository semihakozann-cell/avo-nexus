"""TR: Açık kaynak ek tabloları indir. EN: Download open supplementary tables."""
import urllib.request,concurrent.futures,pathlib,pandas as pd
root=pathlib.Path(__file__).parent
(root/'data').mkdir(exist_ok=True)

def get(n):
 p=root/'data'/f's{n:03d}.xlsx'
 if not p.exists():
  u=f'https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0330523.s{n:03d}&type=supplementary'
  p.write_bytes(urllib.request.urlopen(u,timeout=50).read())
 return p
for p in concurrent.futures.ThreadPoolExecutor(3).map(get,[4,6,8]):
 d=pd.read_excel(p,sheet_name=None)
 print(p.name,p.stat().st_size)
 for k,v in d.items(): print(k,v.shape)
 if p.name=='s008.xlsx': print('COLUMNS',d['Bin_summary'].columns.tolist())
 if p.name=='s004.xlsx':
  for k,v in d.items():
   if 'domain' in k.lower() or k=='Abundances__Bin_ID': print(k,v.head(5).to_string(index=False)); print('SUMS',v.select_dtypes('number').sum().to_dict())
