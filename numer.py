import numpy as np
import matplotlib.pyplot as plt


tp1 = []
tp2 = []
tp3 = []


#a_water = 0.143
#a_water = 0.143
#a_ice = 0.16

a_water = 1.4 * 10 ** -1
a_ice = 1.6 * 10 ** -1
length = 120 #mm
center = (length / 2)
time = 6000 #seconds
nodes = 100
T_0 = 275 # temperature of surrounding water
T_0_st = 275 # temperature of surrounding air
T_pipe = 253 # temperature of pipe
r_pipe = 20 # mm, radius of the pipe
c_p = length / 2 # setting center of pipe

rho = 10 ** -6 # density of water
lmb = 330000 # teplota plavlenia






# Initialization

dx = length / nodes
dy = length / nodes

q = rho * ( dx ** 3) * lmb


rho = 10 ** -6
lmb = 330000
q_phase_change = rho * dx ** 3 * lmb
print("Q PHASE CHANGE IS", q_phase_change)
dt = min(dx**2 / (4 * a_water), dy**2 / (4 * a_water), dx**2 / (4 * a_ice), dy**2 / (4 * a_ice))
#dt = 1
print("DT is", dt)

t_nodes = int(time/dt)

u = np.zeros((nodes, nodes)) + T_0 # Setting T_0
cond_map = np.zeros((nodes, nodes)) # state of mattter (ice or water)

# setting pipe parameters
r_pipe_units = int(r_pipe / dx)
c_p_units = int(c_p / dx) # center of pipe in units (mb we dont need in now)
center_u = center / dx # center of system in units

# setting initial conditions
for x in range (nodes):
    for y in range (nodes):
        if ((x - center_u) ** 2 + (y - center_u) ** 2 > center_u ** 2):
            u[x][y] = T_0_st

def setpipe (r, cx, cy):
    r_units = int(r / dx)
    cx_units = int(cx / dx)
    cy_units = int(cy / dx)
    for x in range(cx_units - r_units, cx_units + r_units):
        for y in range(cy_units - r_units, cy_units + r_units):
            if ((x - cx_units) ** 2 + (y - cy_units) ** 2 <= r_units ** 2):
                u[x][y] = T_pipe
                cond_map[x][y] = 2

def update_cond_map(cond_map):
    cou = 0
    for i in range(1, nodes - 1):
        for j in range(1, nodes - 1):
            if (u[i][j] < 273 and cond_map[i][j] != 2):
                cond_map[i][j] = 1
                cou = cou + 1
    print(cou)


marked_points = []
def boundary(i, j):
    cnd = [0, 0]
    for x in range (i-1, i+2):
        for y in range (j - 1, j + 2):
            if (cond_map[x][y] == 0):
                cnd[0] += 1
            elif (cond_map[x][y] == 1):
                cnd[1] += 1
    if (cnd[0] != cnd[1]):
        print("boundary")
        #plt.scatter(x, y)
        marked_points.append((x, y))
        return True
    else:
        return False

def gettemp(res, x, y):
    tmp = u[int(x / dx), int(y / dx)]
    res.append(tmp)
    print("T: ", tmp)


def read_file(file_path):
    # Initialize an empty array
    array = []

    try:
        # Open the file in read mode
        with open(file_path, 'r') as file:
            # Read each line in the file
            for line in file:
                # Remove any leading/trailing whitespaces and convert the string to an integer
                #number = (line.strip())
                number = float(line)

                # Append the number to the array
                array.append(number)

    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

    return array

# Boundary Conditions

u[0, :] = np.linspace(T_0_st, T_0_st, nodes)
u[-1, :] = np.linspace(T_0_st, T_0_st, nodes)
u[:, 0] = np.linspace(T_0_st, T_0_st, nodes)
u[:, -1] = np.linspace(T_0_st, T_0_st, nodes)


# Visualizing with a plot

fig, axis = plt.subplots()

pcm = axis.pcolormesh(u, cmap=plt.cm.jet, vmin=T_pipe, vmax=T_0)
plt.colorbar(pcm, ax=axis)



# drawing boundary
'''for x in range (nodes):
    for y in range (nodes):
        if ((x - center_u) ** 2 + (y - center_u) ** 2 == center_u ** 2):
            plt.scatter(x, y, color = 'black')
'''

def main():
    counter = 0
    framenum = 0

    while counter < time:


        update_cond_map(cond_map)
        w = u.copy()
        gettemp(tp1,60, 80 )
        gettemp(tp2, 60, 95)
        gettemp(tp3, 60, 110)



        for i in range(1, nodes - 1):
            for j in range(1, nodes - 1):
                setpipe(5, c_p-18 , c_p)
                setpipe(5, c_p+18, c_p)

                t1 = w[i, j]
                t2 = w[i - 1, j]
                t3 = w[i + 1, j]
                t4 = w[i, j - 1]
                t5 = w[i, j + 1]
                T = cond_map[i][j]

                # setting a pipe
                '''setpipe(5, c_p - 18, c_p)
                setpipe(5, c_p + 18, c_p)'''
                dd_ux = (t2 - 2 * t1 + t3) / (dx ** 2)
                dd_uy = (t4 - 2 * t1 + t5) / (dy ** 2)

                r = dt / (dx ** 2)


                if ((t1 < 273 and t2 < 273 and t3 < 273 and t4 < 273 and t5 < 273)):
                    #u[i][j] = r * t1 + ((dt * a_ice) / (dx ** 2)) * (t3 + t2 + t5 + t4 - 4 * t1)
                    u[i][j] = t1 + ((dt * a_ice) / (dx ** 2)) * (t3 + t2 + t5 + t4 -4*t1)
                    marked_points.append((i, j))

                #if (boundary(i, j)):


                else:
                    #u[i, j] = dt * a_water * (dd_ux + dd_uy) + w[i, j]
                    u[i][j] = t1 + ((dt * a_water) / (dx ** 2)) * (t3 + t2 + t5 + t4 -4*t1)






                '''if np.isclose(u[i, j], 273, atol=0.1):
                    plt.scatter(j, i, color='black')
                if (cond_map[i][j] == 1):
                    plt.scatter(j, i, color='black')'''

        counter += dt

        print("t: {:.3f} [s], Average temperature: {:.2f} Kelvins".format(counter, np.average(u)))

        # Updating the plot
        x_marked = [coord[1] for coord in marked_points]
        y_marked = [coord[0] for coord in marked_points]
        #plt.scatter(x_marked, y_marked, color='black', marker='x')


        pcm.set_array(u)
        axis.set_title("Distribution at t: {:.3f} [s].".format(counter))
        fig.savefig(
            f"frames/frame_{framenum}.png"
        )

        framenum += 1
        plt.pause(0.01)


main()

# writing into result files
with open("results_tp1.txt", "w") as f:
    for i in range (len(tp1)):
        f.write(str(tp1[i]))
        f.write("\n")
with open("results_tp2.txt", "w") as f:
    for i in range (len(tp2)):
        f.write(str(tp2[i]))
        f.write("\n")

with open("results_tp3.txt", "w") as f:
    for i in range (len(tp3)):
        f.write(str(tp3[i]))
        f.write("\n")

tp1_f = read_file("results_tp1.txt")
tp2_f = read_file("results_tp2.txt")
tp3_f = read_file("results_tp3.txt")
fig3, ax3 = plt.subplots()

x_tp_1 = np.linspace(0, len(tp1_f), len(tp1_f)) * dt
plt.scatter(x_tp_1, tp1_f)
plt.scatter(x_tp_1, tp2_f)
plt.scatter(x_tp_1, tp3_f)


plt.show()
print("DT : ", dt)