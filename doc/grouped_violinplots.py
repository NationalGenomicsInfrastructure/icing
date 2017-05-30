"""
Grouped violinplots with split violins
======================================

_thumb: .5, .47
"""
import seaborn as sns
import pandas
sns.set(style="whitegrid", palette="pastel", color_codes=True,font_scale=3)
#sns.set(font_scale=3)
#sns.set_context("paper")
# Load the example tips dataset
readLengths = pandas.read_csv("grouped.violin.csv")
#readLengths = pandas.read_csv("grouped.violin.log.csv")

# Draw a nested violinplot and split the violins for easier comparison
ax=sns.violinplot(data=readLengths, split=True,x="barcode",y="readlength", hue="readtype",palette={"ALL":"b","2D":"y"}, saturation=1, scale="count", cut=0)
ax.set_ylim(bottom=0,top=8000)
ax.legend_.remove()
#sns.despine(left=True)
sns.plt.show()
