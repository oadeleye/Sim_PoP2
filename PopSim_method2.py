from __future__ import division
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import time
from matplotlib.pyplot import pause
from math import *
from math import log
import powerlaw


degfile = open('degree.csv', 'w')
matfile = open('matrix.csv', 'w')


def read_nodes():
    nodes = []
    node_file = open('Sorted_services.txt', 'r')
    for node in node_file:
        temp = node.split(';')
        nodes.append(temp[1])
    return (nodes)


# calculate C(n,k)
def comb(n, k):
    return factorial(n) / factorial(k) / factorial(n - k)


# set initial fully connected network as 4 nodes
m0 = 4
init_network_degree = 2 * comb(m0, 2)
pause_time =0.1

#FLAGS
show_degr = 1
draw_powerlaw = True
plot_network = 0


# node position map for plotting
node_pos = {}

# draw the initial network G with random positions
def plot_initial_network(G):
    for i in G.nodes():
       npos = nx.spring_layout(G)

    #Initial plot with only start nodes
    fig = plt.figure('PoP-Sim Network')
    fig.text(0, 0.97, 'starting nodes: green ',style='italic',fontsize=13)
    fig.text(0, 0.94, 'New node: blue ',style='italic',fontsize=13)
    fig.text(0, 0.91, 'Previously added nodes: red', style='italic',fontsize=13)
    nx.draw_networkx(G,npos,node_color = 'green')
    plt.draw()


def plot_new_edges(G, new_edges,i,n):
    plt.clf()
    # Create a color map for nodes to differenctiate between: stating nodes (green), new node (blue) and already added nodes (red)
    color_map = []
    for j in G.nodes():
        if j < m0:
            color_map.append('green')
        elif j == n[m0 + i] :
            color_map.append('blue')
        else: color_map.append('red')
    # Define new node's position and draw the graph
    node_pos= nx.spring_layout(G)
    nx.draw_networkx(G, node_pos, node_color=color_map)
    nx.draw_networkx_edges(G, node_pos,new_edges, width = 2.0 , edge_color = 'b' )
    fig = plt.figure('PoP-Sim Network')
    fig.text(0, 0.97, 'starting nodes: green                 Iteration: '+ str(i+1),style='italic',fontsize=14)
    fig.text(0, 0.94, 'New node: blue ['+str(m0 + i) + ']',style='italic',fontsize=14)
    fig.text(0, 0.91, 'Previously added nodes: red', style='italic',fontsize=14)
    plt.draw()
    pause(pause_time)


def Pop_Sim2( m, gamma):
    n=read_nodes()
    sim_matrix = gen_matrix()
    # initialize graph
    G = nx.Graph()


    # Connect all the initial m0 nodes to the graph
    # need to use the nodelist sorted by birth date
    for i in range(m0):
        G.add_node(n[i])
        for j in G:
            if (n[i] != j):
                G.add_edge(n[i], j)
    # could be removed if plotting is not needed
    if plot_network:
        plot_initial_network(G)

    # Adding new node
    N = len(sim_matrix)

    for i in range(N - m0):
        new_edges = []
        # Compute duration of calculations for consistent timing
        loop_start_time = time.time()
        # select neighbors the new node will connect to according to the  hyperbolic distance
        node_index = choose_neighbour (G,gamma, sim_matrix,m,int(n[i+m0]))
        G.add_node(n[m0 + i])
        # m must not be less than m0
        for node in node_index:
            G.add_edge(n[m0+i],node)
            new_edges.append((n[m0 + i],node))

        degfile.write(str(new_edges) + "\n")

        if plot_network:
            plot_new_edges(G, new_edges, i,n)
    plt.close()

    loop_duration = time.time() - loop_start_time
    # Pause for the needed time, taking the calculation time into account
    if pause_time - loop_duration > 0:
        pause(pause_time - loop_duration)

    if draw_powerlaw:
        gen_degree(G)
    # nx.draw(G)

    if show_degr:
        print('Press any key to continue')
        raw_input()
        powerlaww(G)
        degr(G)
    else:
        print('Press any key to exit')
        raw_input()


def gen_matrix():
    rwr_matrix=np.genfromtxt('real_rwr.csv', delimiter=';')
    np.fill_diagonal(rwr_matrix, 0)  # zero diagonal of rwr matrix
    rwr_matrix=normalize_matrix(rwr_matrix)   # normalize rwr matrix
    uniform_matrix=(1/(len(rwr_matrix)-1))+np.zeros((len(rwr_matrix),len(rwr_matrix))) #generate unifortm matrix
    np.fill_diagonal(uniform_matrix, 0)
    sim_matrix = rwr_matrix + uniform_matrix#  add two matrix
    sim_matrix = (sim_matrix / 2)
    sim_matrix = minus_one(sim_matrix)
    return(sim_matrix)


def normalize_matrix(matrix):
    row_sums = matrix.sum(axis=1)
    new_matrix = matrix / row_sums[:, np.newaxis]
    return new_matrix

def minus_one(matrix):
    l=(len(matrix),len(matrix))
    mat1=np.ones(l)
    np.fill_diagonal(mat1,0)
    new_mat=np.subtract(mat1,matrix)
    return new_mat

def choose_neighbour(G,gamma,sim,m,t):
    neb=new_hyperbolic_dist(G,m,gamma,sim,t)
    print(neb)
    nodes=list(G.nodes())
    T=[nodes[i] for i in neb]
   # print(T)
    return T



def new_hyperbolic_dist(G,m,gamma,sim,t):
    distance=[]
    neigbours=[]
    beta = 1 / (gamma - 1)
    for node in G.nodes():
        rs = (beta * log(int(node)) + (1 - beta) * (log(t)))
        rt=log(t)
        theta =(sim[int(node) - 1][t - 1] * (pi/2))
        d=cosh (rs)* cosh(rt)-sinh(rs) * sinh(rt)*cos(theta)
        if (d<1):
            d=1
            d=acosh(d)
            distance.append(d)
        else:
            d = acosh(d)
            distance.append(d)
    sum_dist = sum(distance)
    distance = [i / sum_dist for i in distance]
    np.asarray(distance)
    np.nan_to_num(distance)
    dist_index = np.argsort(distance)
    for i in range(m):
        neigbours.append(dist_index[i])
    return neigbours

# save degree in a text file for further analysis
def gen_degree(G):
    plt.close()
    degfile = open('degreelist.csv', 'w')
    deg_list = []
    for n in G.nodes():
        deg_list.append(G.degree(n))
    degrees = np.array(deg_list)
    for d in deg_list:
        degfile.write(str(d) + '\n')
    return (degrees)


# powerlaw fitting and degree exponent.
def powerlaww(G):
    degree_list = gen_degree(G)
    fit = powerlaw.Fit(degree_list, linear_bins=False)
    print("The degree of exponent is:", fit.power_law.alpha)
    print("gamma's statndard error is:", fit.power_law.sigma)
    print("Likelihood of PL and Exponetial:", fit.distribution_compare('power_law', 'exponential'))
    fig1 = fit.plot_pdf(color='r', linewidth=0, marker='o')
    #fit.plot_ccdf(color='b', linewidth=2)
    #fit.power_law.plot_pdf(color='b', linestyle='--')
    # fit.plot_ccdf(color='r', linewidth=2,ax=fig2)
    #   fit.power_law.plot_ccdf(color='r',linestyle='--',ax=fig2)
    #   powerlaw.plot_pdf(rands,linear_bins=True,color='r')
    plt.pause(100)


def degr(G, scale='log',alpha=.8,low=1,high=1,ec=1):
    plt.close()
    num_nodes=G.number_of_nodes()
    max_degree=0
    for n in G.nodes():
        if G.degree(n)>max_degree:
            max_degree=G.degree(n)
    x=[]
    y_tmp=[]
    for i in range (max_degree+1):
        x.append(i)
        y_tmp.append(0)
        for n in G.nodes():
            if G.degree(n)==i:
                y_tmp[i] += 1
    y=[i/num_nodes for i in y_tmp]
    deg, =plt.plot(x,y,label='Degree Distribution',linewidth=0 , marker='o', markersize=7, color='blue',alpha=alpha)
    # Check for the lin / log parameter and set axes scale
    if scale == 'log':
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Degree distribution (log-log scale)')
        # add theoretical distribution line k^-3
        w = [a for a in range(low,high)]
        z = []
        for i in w:
            x = (i ** -3) * ec  # set line's length and fit intercept
            z.append(x)
        plt.plot(w, z, 'k-', color='#7f7f7f')
    else:
        plt.title('Degree distribution (linear scale)')

    plt.ylabel('Pk')
    plt.xlabel('k')
    plt.title( 'Degree Distribution (Closest Node Method )')
    plt.show()


if __name__ == "__main__":

    Pop_Sim2(2, 3) # inputs are m and gamma respectively3
