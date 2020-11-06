import matplotlib.pyplot as plt

x = 1
y = 2

plt.figure(1)
# plt.plot(x,y, label='Parte Superior')
plt.plot(x, y, marker='x', color='r', label='CLmax')
plt.legend()
plt.show()