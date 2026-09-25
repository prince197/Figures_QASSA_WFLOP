import numpy as np, truth_eval as te, json, sys
d,o=sys.argv[1],sys.argv[2]
t=np.load(f'{d}/truth.npz'); c=te.load_curve(f'{d}/env/power_curve.csv')
xy=np.loadtxt(f'{o}/layout.csv',delimiter=',',skiprows=1)[:,1:]
P=te.farm_power(xy,t['U'],t['theta'],c); Pf=te.free_power(t['U'],c)
y=json.load(open(f'{o}/yield.json'))
print(json.dumps(dict(true_net=P,eff=P/len(xy)/Pf,pred_err_pct=100*(y['net_mean_power_kw']/P-1),gross_err_pct=100*(y['gross_mean_power_kw']/(len(xy)*Pf)-1))))
