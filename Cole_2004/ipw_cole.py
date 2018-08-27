########################################################################################################################
# Demonstration of zEpid Python 3 library
#       Inverse Probability Weights) Demonstration of the IPTW class, using an example from Cole and Hernan. "Adjusted
#                   survival curves with inverse probability weights"
#                   DOI: 10.1016/j.cmpb.2003.10.004
#
# Using this code:
#   Run code through the terminal. The Ewing Sarcoma data will be loaded from zEpid and the process of estimating IPW
#       will be completed
#
#   Graph of survival curves (weighted and unweighted) are created
#
# Versions:
#   Python 3.6.3
#   numpy v1.14.5
#   pandas v0.23.0
#   lifelines v0.14.6
#
########################################################################################################################
# Code by: Paul Zivich                                                      Last edit: 2018/8/7
########################################################################################################################

import pandas as pd
from lifelines import KaplanMeierFitter
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import zepid as ze
from zepid.causal.ipw import IPTW

df = ze.load_ewing_sarcoma_data()
df['id'] = [i for i in range(df.shape[0])]
ipt = IPTW(df,  # pandas dataframe object
           treatment='treat',  # treatment/exposure column label
           stabilized=True)  # optional argument with True as default
ipt.regression_models('ldh')
ipt.fit()  # generating IPTW (will be contained in IPTW.Weight)
df['iptw'] = ipt.Weight  # extracting weights and adding to a new column in our dataframe
print(round(df['iptw'].sum(), 5))  # Checking that stabilized weights add to original sample size

df_t = df.loc[df['treat'] == 1].copy()  # Dividing data into novel treatment vs standard treatment
df_u = df.loc[df['treat'] == 0].copy()

# Unweighted Kaplan-Meier
kmn_t = KaplanMeierFitter()  # Among Novel treatment
kmn_t.fit(durations=df_t['time'], event_observed=df_t['outcome'])
unad_t = kmn_t.event_table
unad_t['survival'] = kmn_t.survival_function_
plot_u_t = unad_t.loc[unad_t['observed'] == 1].copy()

kmn_u = KaplanMeierFitter()  # Among Standard treatment
kmn_u.fit(durations=df_u['time'], event_observed=df_u['outcome'])
unad_u = kmn_u.event_table
unad_u['survival'] = kmn_u.survival_function_
plot_u_u = unad_u.loc[unad_u.index < 1500].copy()

# Adjusted Kaplan-Meier
kmn_t = KaplanMeierFitter()  # Among Novel treatment
kmn_t.fit(durations=df_t['time'], event_observed=df_t['outcome'], weights=df_t['iptw'])
ad_t = kmn_t.event_table
ad_t['survival'] = kmn_t.survival_function_
plot_a_t = ad_t.loc[ad_t['observed'] != 0].copy()

kmn_u = KaplanMeierFitter()  # Among Standard treatment
kmn_u.fit(durations=df_u['time'], event_observed=df_u['outcome'], weights=df_u['iptw'])
ad_u = kmn_u.event_table
ad_u['survival'] = kmn_u.survival_function_
plot_a_u = ad_u.loc[ad_u.index < 1500].copy()

# Recreating Figure 1
plt.figure(figsize=[8, 4])
gspec = gridspec.GridSpec(1, 2)
una = plt.subplot(gspec[0, 0:1])
adj = plt.subplot(gspec[0, 1:])

una.step(plot_u_t.index, plot_u_t.survival, linestyle='--', c='gray', where='post', label='Novel Tx')
una.step(plot_u_u.index, plot_u_u.survival, c='k', where='post', label='Standard Tx')
una.set_ylim([0, 1.01])
una.set_xlim([0, 1500])
una.set_yticks([0.0, 0.25, 0.50, 0.75, 1.00])
una.set_xticks([0, 500, 1000, 1500])
una.legend()
una.set_title('Unadjusted')
adj.step(plot_a_t.index, plot_a_t.survival, linestyle='--', c='gray', where='post', label='Novel Tx')
adj.step(plot_a_u.index, plot_a_u.survival, c='k', where='post', label='Standard Tx')
adj.set_ylim([0, 1.01])
adj.set_xlim([0, 1500])
adj.set_yticks([0.0, 0.25, 0.50, 0.75, 1.00])
adj.set_xticks([0, 500, 1000, 1500])
adj.legend()
adj.set_title('Weighted')
plt.show()

###############################################################################
# Calculating Inverse Probability of Censoring Weights
from zepid.causal.ipw import IPCW

ipc = IPCW(df,  # pandas dataframe object
           idvar='id',  # unique ID for participant
           time='time',  # time column label
           event='outcome',  # outcome column label
           flat_df=True)  # whether a flat df is input. Converts to a long dataframe behind the scenes
model_n = 't_enter + I(t_enter**2) + I(t_enter**3)'
model_d = 'treat + ldh + t_enter + + I(t_enter**2) + I(t_enter**3)'
ipc.regression_models(model_denominator=model_d,  # dependent model form for denominator
                      model_numerator=model_n)  # dependent model form for numerator
ipc.fit()  # estimating IPCW
df['ipcw'] = ipc.Weight  # extracting weights and adding to a new column in our dataframe
# TODO still need to correct this bit

ipc.df.info()

# Repeating but with IPTW*IPCW
df['w'] = df['iptw'] * df['ipcw']
df_t = df.loc[df['treat'] == 1].copy()  # Dividing data into novel treatment vs standard treatment
df_u = df.loc[df['treat'] == 0].copy()

kmn_t = KaplanMeierFitter()  # Among Novel treatment
kmn_t.fit(durations=df_t['time'], event_observed=df_t['outcome'], weights=df_t['w'])
ad_t = kmn_t.event_table
ad_t['survival'] = kmn_t.survival_function_
plot_a_t = ad_t.loc[ad_t['observed'] != 0].copy()

kmn_u = KaplanMeierFitter()  # Among Standard treatment
kmn_u.fit(durations=df_u['time'], event_observed=df_u['outcome'], weights=df_u['w'])
ad_u = kmn_u.event_table
ad_u['survival'] = kmn_u.survival_function_
plot_a_u = ad_u.loc[ad_u.index < 1500].copy()

# Recreating Figure 1
plt.figure(figsize=[8, 4])
gspec = gridspec.GridSpec(1, 2)
una = plt.subplot(gspec[0, 0:1])
adj = plt.subplot(gspec[0, 1:])

una.step(plot_u_t.index, plot_u_t.survival, linestyle='--', c='gray', where='post', label='Novel Tx')
una.step(plot_u_u.index, plot_u_u.survival, c='k', where='post', label='Standard Tx')
una.set_ylim([0, 1.01])
una.set_xlim([0, 1500])
una.set_yticks([0.0, 0.25, 0.50, 0.75, 1.00])
una.set_xticks([0, 500, 1000, 1500])
una.legend()
una.set_title('Unadjusted')
adj.step(plot_a_t.index, plot_a_t.survival, linestyle='--', c='gray', where='post', label='Novel Tx')
adj.step(plot_a_u.index, plot_a_u.survival, c='k', where='post', label='Standard Tx')
adj.set_ylim([0, 1.01])
adj.set_xlim([0, 1500])
adj.set_yticks([0.0, 0.25, 0.50, 0.75, 1.00])
adj.set_xticks([0, 500, 1000, 1500])
adj.legend()
adj.set_title('Weighted')
plt.show()
