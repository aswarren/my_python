"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""
import numpy as np
import matplotlib.pyplot as plt


n_groups = 5

means_google = (20, 35, 30, 35, 27)
std_google = (2, 3, 4, 1, 2)

means_actual = (25, 32, 34, 20, 25)
std_actual = (3, 5, 2, 3, 3)

fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.35

opacity = 0.4
error_config = {'ecolor': '0.3'}

rects1 = plt.bar(index, means_google, bar_width,
                 alpha=opacity,
                 color='g',
                 yerr=std_google,
                 error_kw=error_config,
                 label='Google')

rects2 = plt.bar(index + bar_width, means_actual, bar_width,
                 alpha=opacity,
                 color='b',
                 yerr=std_actual,
                 error_kw=error_config,
                 label='Actual')

plt.xlabel('Group')
plt.ylabel('Scores')
plt.title('Scores by group and gender')
plt.xticks(index + bar_width, ('Fifth Element', 'Proteomics', 'Genomics', 'Lab Equipment', 'Headphones'))
plt.legend()

plt.tight_layout()
plt.show()
