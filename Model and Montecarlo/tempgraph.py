import matplotlib.pyplot as plt
legend_drawn_flag = True

g0="Pd3"
g1="Au1Pd2"
g2="Au2Pd1"
g3="Au3"

SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

# line 1 points
x1 = [300, 400, 500,600,700,800,900,1000,1100,1200]
y1 = [0.00049799520226218, 0.001956123621530463, 0.003699900419009635,0.006592597506174551,0.00933079394497886,0.010239516621553976,0.0141706429579085,0.017857704278338198,0.019363903210177322,0.02147071778830238]
# plotting the line 1 points
plt.plot(x1, y1, label=g0.translate(SUB),linewidth=3,color='purple')

# line 2 points
x2 = [300, 400, 500,600,700,800,900,1000,1100,1200]
y2 = [0.08930553779881538, 0.11389797792100619, 0.13503785672597043,0.15433546043899077,0.173675770139508,0.18473600192433776,0.20424718522116525,0.217520891988594,0.2239064677546195,0.2379630285229232]
# plotting the line 2 points
plt.plot(x2, y2, label=g1.translate(SUB),linewidth=3,color='red')


x3 = [300, 400, 500,600,700,800,900,1000,1100,1200]
y3 = [0.5828296596199959,0.5570384574901752,0.5556586615017649,0.5486624532898293,0.5380778111347297,0.5335432698050266,0.5298238369902671,0.5203176677353335,0.5092321538882405,0.5079345777446007]
# plotting the line 2 points
plt.plot(x3, y3, label=g2.translate(SUB),linewidth=3,color='blue')


x4 = [300, 400, 500,600,700,800,900,1000,1100,1200]
y4 = [0.3233323439400727,0.32578530585056303,0.3095703856906148,0.2926392753261108,0.28174499492136734,0.2697443152849667,0.24792371808539446,0.24421231683553396,0.2353160799155334,0.2253958293879351]
# plotting the line 2 points
plt.plot(x4, y4, label=g3.translate(SUB),linewidth=3,color='green')

# naming the x axis
plt.xlabel('Temperature (K)',fontsize=24)
# naming the y axis
plt.ylabel('Active Site Distribution (%)',fontsize=24)
# giving a title to my graph
# plt.title('Distribution of Active sites across Temperatures')

# show a legend on the plot
plt.legend(loc=0, frameon=legend_drawn_flag,fontsize=24)
plt.xticks(fontsize=24)
plt.yticks(fontsize=24)

# function to show the plot
plt.show()






# print('space')
# plt.plot(gl0)
# plt.plot(gl1)
# plt.plot(gl2)
# plt.plot(gl3)
# legend_drawn_flag = True
# plt.legend(["g0", "g1",'g2','g3'], loc=0, frameon=legend_drawn_flag)
# plt.show()