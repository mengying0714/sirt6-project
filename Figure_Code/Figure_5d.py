###############################################################################
# Code for Figure 5d
# Language: Python 3
# Required package: matplotlib
#Mengying Zhang (zhangmyharper@gmail.com) 
###############################################################################

import matplotlib.pyplot as plt

counts = [13900, 1284]
categories = ["Old-Gained\n(SIRT6 reversed)", "Old-Depleted\n(SIRT6 restored)"]
colors = ["#E05C5C", "#5C8AE0"]

break_low = 1550
break_high = 5000

fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(4.5, 6), sharex=True, gridspec_kw={"height_ratios": [1.8, 1.5], "hspace": 0.05})

for i, (cnt, col) in enumerate(zip(counts, colors)):
    ax_top.bar(i, cnt, color=col, width=0.5, zorder=3)

ax_top.set_ylim(break_high, 14250)
ax_top.set_yticks([5000, 8000, 11000, 14000])
ax_top.set_yticklabels(["5,000", "8,000", "11,000", "14,000"], fontsize=15)
ax_top.spines["top"].set_visible(False)
ax_top.spines["bottom"].set_visible(False)
ax_top.spines["right"].set_visible(False)
ax_top.tick_params(axis="x", bottom=False, labelbottom=False)

for i, (cnt, col) in enumerate(zip(counts, colors)):
    ax_bot.bar(i, cnt, color=col, width=0.5, zorder=3)

ax_bot.set_ylim(0, break_low)
ax_bot.set_yticks([0, 500, 1000, 1500])
ax_bot.set_yticklabels(["0", "500", "1,000", "1,500"], fontsize=15)
ax_bot.spines["top"].set_visible(False)
ax_bot.spines["right"].set_visible(False)

d = 0.015
ax_top.plot([-d, d], [-d, d], transform=ax_top.transAxes, color="black", linewidth=1.2, clip_on=False)
ax_bot.plot([-d, d], [1-d, 1+d], transform=ax_bot.transAxes, color="black", linewidth=1.2, clip_on=False)

ax_bot.set_xticks([0, 1])
ax_bot.set_xticklabels(categories, fontsize=15, fontweight="normal")

fig.text(0.02, 0.5, "Number of chromatin loops", va="center", rotation="vertical", fontsize=15, fontweight="normal")
fig.subplots_adjust(left=0.25)

plt.savefig("loop_category_barplot_break5000.pdf", dpi=300, bbox_inches="tight")
plt.savefig("loop_category_barplot_break5000.png", dpi=300, bbox_inches="tight")
plt.show()