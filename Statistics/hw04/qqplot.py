'''
QQ plots

'''


def qqplots(res):
  from scipy import stats
  import matplotlib.pyplot as plt
  import numpy as np

  res.sort()
  n = len(res)
  p = np.linspace(0 + 1./(n-1), 1-1./(n-1), n)
  quants = np.zeros_like(res)
  for i in range(n):
    quants[i] = stats.scoreatpercentile(res, p[i]*100)

  mu = res.mean()
  sigma = res.std()
  y = stats.norm.ppf(p, loc=0, scale=1)

  figx = plt.figure(figsize=(10,6))

  plt.scatter(y, quants)


  plt.title('Normal - Quantile Plot')
  plt.ylabel('Anual rainfall')
  plt.xlabel('Quantiles of N(0,1)')

