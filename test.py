import matplotlib.pyplot as plt
# doubling the width of markers
x = [0,2,4,6,8,10,12,14,16,18,20,22,24]
y = [0]*len(x)
s = [400*n for n in range(len(x))]
plt.scatter(x,y,s=s)
plt.show()