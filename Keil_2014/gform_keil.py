########################################################################################################################
# Demonstration of zEpid Python 3 library
#       Time-varying g-formula) Demonstration of the TimeVaryGFormula class, using an example from Keil et al. "The
#                   parametric G-formula for time-to-event data: towards intuition with a worked example"
#                   DOI: 10.1097/EDE.0000000000000160
#
# Using this code:
#   Run code through the terminal. The GvHD will be loaded from zEpid and the process of fitting the model
#       will be completed. Note, this code may take awhile unless the resamples are maximum time are
#       changed.
#
#   Two graphs (observed vs natural & natural vs prevent GvHD) are created. The Cox PH model results are printed to the
#       the terminal
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

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
import matplotlib.pyplot as plt
import zepid as ze
from zepid.causal.gformula import TimeVaryGFormula

df = ze.load_gvhd_data()
df['tomorrow'] = df['day']+1

# G-Compuation Algorithm Process
g = TimeVaryGFormula(df,  # pandas dataframe object containing data to fit to
                     idvar='id',  # unique ID for participants
                     exposure='gvhd',  # exposure column label (gvhd)
                     outcome='d',  # outcome column label (d)
                     time_in='day',  # entry time column label
                     time_out='tomorrow')  # exit time column label

# Fitting exposure model to observed data
exp_m = ('all + cmv + male + age + agecurs1 + agecurs2 + platnormm1 + daysnoplatnorm + relapsem1 + daysnorelapse + '
         'day + daysq + wait')
g.exposure_model(exp_m,  # dependent variable model to predict exposure
                 restriction="g['gvhdm1']==0")  # restricting for ITT

# Fitting outcome model to observed data
out_m = ('all + cmv + male + age  + agesq + day + daysq + gvhd + platnorm + daysnoplatnorm + relapse + daysnorelapse +'
         ' daycu + wait + day*gvhd + daysq*gvhd + daycu*gvhd')
g.outcome_model(out_m)  # dependent variable model to predict outcome

# Fitting model for PlatNorm
prp_m = 'all + cmv + male + age + agecurs1 + agecurs2 + gvhdm1 + daysgvhd + daysnorelapse + wait'
prp_recode = ("g['platnorm'] = np.where(g['platnormm1']==1,1,g['platnorm']);"
              "g['daysnoplatnorm'] = np.where(g['platnorm']==0,g['daysnoplatnorm']+1,g['daysnoplatnorm']);"
              "g['daysplatnorm'] = np.where(g['platnorm']==0,g['daysplatnorm'],g['daysplatnorm']+1);")
g.add_covariate_model(label=1,  # order to fit the covariate model
                      covariate='platnorm',  # covariate to predict
                      model=prp_m,  # dependent variable model to predict covariate (platnorm)
                      recode=prp_recode,  # recoding process to execute each loop
                      var_type='binary',  # covariate variable type
                      restriction="g['platnormm1']==0")  # restricting for ITT

# Fitting model for probability of relapse
prr_m = ('all + cmv + male + age + agecurs1 + agecurs2 + gvhdm1 +  daysgvhd + platnormm1 + daysnoplatnorm + day + '
         'daysq + wait')
relapse_recode = ("g['relapse'] = np.where(g['relapsem1']==1,1,g['relapse']);"
                  "g['daysnorelapse'] = np.where(g['relapse']==0,g['daysnorelapse']+1,g['daysnorelapse']);"
                  "g['daysrelapse'] = np.where(g['relapse']==0,g['daysrelapse'],g['daysrelapse']+1);")
g.add_covariate_model(label=2,  # order to fit the covariate model
                      covariate='relapse',  # covariate to predict
                      model=prr_m,  # dependent variable model to predict covariate (relapse)
                      recode=relapse_recode,  # recoding process to execute each loop
                      var_type='binary',  # covariate variable type
                      restriction="g['relapsem1']==0")  # restricting for ITT

# Fitting model for censoring
cen_m = 'all + cmv + male + age + agesq + daysgvhd + daysnoplatnorm + daysnorelapse + day + daysq + daycu + wait'
g.add_covariate_model(label=3,  # order to fit the covariate model
                      covariate='censlost',  # covariate to predict
                      model=cen_m,  # dependent variable model to predict covariate (relapse)
                      var_type='binary')  # covariate variable type

# Fitting the g-formula
g.fit(treatment="((g['gvhd']==1) | (g['gvhdm1']==1))",  # treatment strategy (ITT natural course)
      lags={'platnorm': 'platnormm1',  # lagging certain variables at end of each loop
            'relapse': 'relapsem1',
            'gvhd': 'gvhdm1'},
      sample=13700,  # number to resample from the population for Monte Carlo
      t_max=1825,  # maximum time to simulate until (5 years)
      in_recode=("g = g.loc[g['censlost']==0].copy();"  # variable recoding to execute at start of loop
                 "g['daysq'] = g['day']**2;"
                 "g['daycu'] = g['day']**3;"),
      out_recode=("g['daysnogvhd'] = np.where(g['gvhd']==0,g['daysnogvhd']+1,g['daysnogvhd']);"  # recoding at end
                  "g['daysgvhd'] = np.where(g['gvhd']==0,g['daysgvhd'],g['daysgvhd']+1);"))

gf = g.predicted_outcomes
gf['d'] = np.where(gf['censlost'] == 1, 0, gf['d'])
gfs = gf.loc[gf.uid_g_zepid != gf.uid_g_zepid.shift(-1)].copy()
kmn = KaplanMeierFitter()
kmn.fit(durations=gfs['day'], event_observed=gfs['d'])

rfn = pd.DataFrame()  # Creating df for CoxPHFitter
rfn['event'] = gfs['d']
rfn['time'] = gfs['day']
rfn['exp'] = 0

g.fit(treatment="none",
      lags={'platnorm': 'platnormm1',
            'relapse': 'relapsem1',
            'gvhd': 'gvhdm1'},
      sample=13700,  # 137000 in original paper, reduced it by a factor of 10
      t_max=1825,  # 1825 in original paper, keeping same since 5 years out
      in_recode=("g['daysq'] = g['day']**2;"
                 "g['daycu'] = g['day']**3;"),
      out_recode=("g['daysnogvhd'] = np.where(g['gvhd']==0,g['daysnogvhd']+1,g['daysnogvhd']);"
                  "g['daysgvhd'] = np.where(g['gvhd']==0,g['daysgvhd'],g['daysgvhd']+1);")) # DONE

gf = g.predicted_outcomes.loc[g.predicted_outcomes['censlost'] != 1].copy()
gfs = gf.loc[gf.uid_g_zepid != gf.uid_g_zepid.shift(-1)].copy()
kmp = KaplanMeierFitter()
kmp.fit(durations=gfs['day'], event_observed=gfs['d'])

rfp = pd.DataFrame()  # Creating df for CoxPHFitter
rfp['event'] = gfs['d']
rfp['time'] = gfs['day']
rfp['exp'] = 1
rf = rfn.append(rfp, ignore_index=True, sort=False)

df.sort_values(by=['id', 'day'])
dfs = df.loc[df.id != df.id.shift(-1)].copy()
kmo = KaplanMeierFitter()
kmo.fit(durations=dfs['day'], event_observed=dfs['d'])

plt.step(kmo.event_table.index, kmo.survival_function_, c='gray', where='post', label='Observed')
plt.step(kmn.event_table.index, kmn.survival_function_, c='firebrick', where='post', label='Natural')
plt.step(kmp.event_table.index, kmp.survival_function_, c='blue', where='post', label='Prevent GvHD')
plt.legend()
plt.show()

cph = CoxPHFitter()
cph.fit(rf, duration_col='time', event_col='event')
cph.print_summary()
print(cph.hazards_)