"""TR: Tuz Gölü kalite eşiği duyarlılık pilotu. EN: Tuz Lake quality-threshold sensitivity pilot.
Run / Çalıştır: python analyze.py
Inputs / Girdiler: PLOS ONE 10.1371/journal.pone.0330523, S1/S3/S5.
Not RNA-seq; no causal inference / RNA-seq değildir; nedensel çıkarım yapılmaz.
"""
from pathlib import Path
import json,hashlib
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=Path(__file__).parent
out=root/'results';out.mkdir(exist_ok=True)
a=pd.read_excel(root/'data/s004.xlsx',sheet_name='Abundances__Bin_ID')
q=pd.read_excel(root/'data/s008.xlsx',sheet_name='Bin_summary')
sm=pd.read_excel(root/'data/s006.xlsx',sheet_name='SM_abundances_taxonomy_phylum')
assert a.Bin_ID.is_unique and q.Bin_ID.is_unique
m=a.merge(q,on='Bin_ID',how='left',validate='one_to_one',indicator=True)
assert m['_merge'].eq('both').all()
assert m[['Completeness','Contamination','Tuz Gölü']].notna().all().all()
assert np.isclose(m['Tuz Gölü'].sum(),100)
assert (m['Tuz Gölü']>=0).all()
assert m['Is_representative'].eq('yes').all()
m['domain_group']=m.Domain.map({'d__Archaea':'Archaea','Unclassified Archaea':'Archaea','d__Bacteria':'Bacteria','Unclassified Bacteria':'Bacteria','Unclassified':'Unclassified'})
assert m.domain_group.notna().all()
w=m['Tuz Gölü']
# Labels describe thresholds, not formal MIMAG quality classes.
# Etiketler eşikleri tanımlar; resmî MIMAG kalite sınıfı değildir.
filters={
 'All representative bins':pd.Series(True,index=m.index),
 'Completeness >=50, contamination <=10':(m.Completeness>=50)&(m.Contamination<=10),
 'Completeness >=70, contamination <=5':(m.Completeness>=70)&(m.Contamination<=5),
 'Completeness >=90, contamination <=5':(m.Completeness>=90)&(m.Contamination<=5),
 'Paper score: completeness -3*contamination >=50':m['Quality_completeness_minus_3xcontamination']>=50,
}
rows=[];profiles=[]
for label,mask in filters.items():
 total=w[mask].sum();arch=w[mask&m.domain_group.eq('Archaea')].sum();bact=w[mask&m.domain_group.eq('Bacteria')].sum()
 rows.append({'filter':label,'passing_bins_all_sites':int(mask.sum()),'passing_bins_Tuz_nonzero':int((mask&(w>0)).sum()),'passing_bins_Tuz_ge_0.1pct':int((mask&(w>=0.1)).sum()),'retained_pct_of_original_bin_abundance':float(total),'excluded_pct':float(100-total),'archaea_pct_of_original':float(arch),'bacteria_pct_of_original':float(bact),'other_unclassified_pct_of_original':float(total-arch-bact),'archaea_pct_within_retained':float(100*arch/total) if total else None})
 p=m.loc[mask].groupby('Phylum')['Tuz Gölü'].sum().sort_values(ascending=False)
 profiles.extend({'filter':label,'phylum':k,'pct_of_original_bin_abundance':float(v),'pct_within_retained':float(v/total*100)} for k,v in p.items())
s=pd.DataFrame(rows);s.to_csv(out/'quality_sensitivity.csv',index=False)
pd.DataFrame(profiles).to_csv(out/'phylum_profiles.csv',index=False)
m.drop(columns='_merge').to_csv(out/'joined_bins.csv',index=False)
# SingleM is normalized separately; denominators differ from bin-derived abundances.
# SingleM ayrı normalize edilir; genom-kümesi bolluklarıyla paydaları farklıdır.
sm['domain']=sm.iloc[:,0].str.split(';').str[0]
sm_sum=float(sm['Tuz Gölü'].sum());domain_sm=(sm.groupby('domain')['Tuz Gölü'].sum()/sm_sum*100).to_dict()
comparison=[]
base=m.groupby('Phylum')['Tuz Gölü'].sum()
sm_p=sm.set_index(sm.columns[0])['Tuz Gölü'];sm_p.index=sm_p.index.str.split(';').str[-1]
sm_p=sm_p.groupby(level=0).sum()/sm_sum*100
for k in sorted(set(base.index)|set(sm_p.index)):
 comparison.append({'phylum':k,'bin_method_relative_pct':float(base.get(k,0)),'SingleM_relative_pct':float(sm_p.get(k,0))})
c=pd.DataFrame(comparison);c['difference_percentage_points']=c.bin_method_relative_pct-c.SingleM_relative_pct;c.to_csv(out/'methods_phylum_comparison.csv',index=False)
hq=pd.read_excel(root/'data/s004.xlsx',sheet_name='HQ_Abundances__Bin_ID')
paper_mask=filters[list(filters)[-1]]
assert set(m.loc[paper_mask,'Bin_ID'])==set(hq.Bin_ID)
# Cross-check direct weighted calculations with aggregated phylum values.
for label,mask in filters.items():
 assert np.isclose(sum(x['pct_of_original_bin_abundance'] for x in profiles if x['filter']==label),w[mask].sum())
labels=['Tüm kümeler','Bütünlük ≥50 / kont. ≤10','Bütünlük ≥70 / kont. ≤5','Bütünlük ≥90 / kont. ≤5','Makalenin kalite puanı ≥50']
fig,axes=plt.subplots(1,2,figsize=(13,5.5),gridspec_kw={'width_ratios':[1.4,1]})
colors=['#178b89','#d4a44b','#94a3b8','#e9edf2']
left=np.zeros(len(s))
for col,label,color in zip(['archaea_pct_of_original','bacteria_pct_of_original','other_unclassified_pct_of_original','excluded_pct'],['Arkea / Archaea','Bakteri / Bacteria','Diğer / sınıflanmamış','Filtreyle dışarıda kalan'],colors):
 vals=s[col].values;axes[0].barh(labels,vals,left=left,label=label,color=color);left+=vals
axes[0].invert_yaxis();axes[0].set_xlim(0,100);axes[0].set_xlabel('Başlangıç genom-kümesi bolluğunun yüzdesi (%)');axes[0].set_title('Kalite filtresi hangi bolluğu dışarıda bırakıyor?')
x=np.arange(2);a0=s.iloc[0];vals_a=[a0.archaea_pct_of_original,domain_sm.get('d__Archaea',0)];vals_b=[a0.bacteria_pct_of_original,domain_sm.get('d__Bacteria',0)]
axes[1].bar(x,vals_a,color=colors[0],label='Arkea');axes[1].bar(x,vals_b,bottom=vals_a,color=colors[1],label='Bakteri');axes[1].bar(x,100-np.array(vals_a)-np.array(vals_b),bottom=np.array(vals_a)+np.array(vals_b),color=colors[2]);axes[1].set_xticks(x,['Genom kümeleri','SingleM']);axes[1].set_ylim(0,100);axes[1].set_ylabel('Her yöntemin kendi içindeki yüzde (%)');axes[1].set_title('Yöntemler farklı profiller verebilir')
for i,v in enumerate(vals_a):axes[1].text(i,v/2,f'{v:.1f}% arkea',ha='center',color='white',fontweight='bold')
fig.suptitle('AVO pilotu: Tuz Gölü / Tuz Lake',fontsize=16,fontweight='bold');fig.text(.5,.015,'DNA metagenomiği • Payda: tüm çevresel DNA değil • Tek saha profili; nedensel/istatistiksel grup karşılaştırması yok',ha='center',fontsize=9)
fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.065),ncol=4,fontsize=9);fig.tight_layout(rect=[0,.18,1,.92]);fig.savefig(out/'tuz_golu_comparison.png',dpi=180,bbox_inches='tight');plt.close(fig)
meta={'source_doi':'10.1371/journal.pone.0330523','table_files':{p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in (root/'data').glob('*.xlsx')},'joins':{'abundance_rows':len(a),'quality_rows':len(q),'matched_abundance_rows':len(m),'quality_rows_not_in_abundance':len(q)-len(m)},'SingleM_input_sum':sm_sum,'SingleM_domain_percent':domain_sm,'source_HQ_bin_set_matches_score_filter':True,'thresholds_are_exploratory':True,'domain_mapping_note':'Unclassified Archaea and Unclassified Bacteria retain their explicitly supplied broad domain. Unassigned phylum remains unassigned.'}
(out/'provenance.json').write_text(json.dumps(meta,indent=2,ensure_ascii=False))
print(s.round(4).to_string(index=False));print('SINGLEM',domain_sm,'sum',sm_sum);print('TOP_METHOD_DIFFERENCES');print(c.reindex(c.difference_percentage_points.abs().sort_values(ascending=False).index).head(8).round(3).to_string(index=False));print('QA',meta['joins'])
