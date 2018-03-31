from __future__ import division
import numpy as np
from math import factorial
import networkx as nx
import csv
import matplotlib.pyplot as plt
import time
from matplotlib.pyplot import pause
from math import *
from math import log
import pandas as pd
import powerlaw


degfile=open('degree.csv','w')
matfile=open('matrix.csv','w')
def read_nodes():
    nodes=[]
    node_file= open('Sorted_services.txt','r')
    for node in node_file:
        temp=node.split(';')
        nodes.append(temp[1])
    return (nodes)

# calculate C(n,k)
def comb(n,k):
    return factorial(n) / factorial(k) / factorial(n - k)

        
# set initial fully connected network as 4 nodes
m0 = 4
init_network_degree = 2*comb(m0,2)
pause_time = 0.1

# Flags - {show degree, plot powerlaw, plot network respectively}
show_degr = 1
draw_powerlaw = 1
plot_network=0

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
    nx.draw_networkx(G,npos, node_color = 'green')
    plt.draw()

def plot_new_edges(G, new_edges, i,n):
    plt.clf()
    # Create a color map for nodes to differenctiate between: stating nodes (green), new node (blue) and already added nodes (red)
    color_map = []
    for j in G.nodes():
        if int(j) < m0:
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

def Pop_Sim(m,gamma):
    sim_matrix=gen_matrix()
    n = read_nodes() # read API nodes
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

    # Adding new nodes
    N=len(sim_matrix)
    for i in range(N-m0):
        # Compute duration of calculations for consistent timing
        loop_start_time = time.time()
       # select neighbors the new node will connect to according to the  hyperbolic distance
        neighbors = choose_neighbour(G,gamma,m,sim_matrix,int(n[i+m0]))
        # A Check to make sure the correct number of neighbors are chosen
        if (len(neighbors) != m):
            print ("Error, number of neighbors is not as expected")
            return
        # Add the new node to the graph
        G.add_node(n[m0 + i])
        # Save new edges in a list for drawing purposed
        new_edges = []

        for nb in neighbors:
            G.add_edge(n[m0 + i], nb)
            new_edges.append((n[m0 + i], nb))
        degfile.write(str(new_edges)+'\n')


        if plot_network:
            plot_new_edges(G, new_edges, i,n)

    plt.close()

    loop_duration = time.time() - loop_start_time
    # Pause for the needed time, taking the calculation time into account
    if pause_time - loop_duration > 0:
        pause(pause_time - loop_duration)

    if draw_powerlaw:
        gen_degree(G)
    #nx.draw(G)

    if show_degr:
        print ('Press any key to continue')
        raw_input()
        powerlaww(G)
        degr(G)

    else:
        print ('Press any key to exit')
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
    normalize_matrix(sim_matrix)
    return(sim_matrix)

def normalize_matrix(matrix):
    row_sums = matrix.sum(axis=1)
    new_matrix = matrix/row_sums[:, np.newaxis]
    return(new_matrix)

def minus_one(matrix):
    l=(len(matrix),len(matrix))
    mat1=np.ones(l)
    np.fill_diagonal(mat1,0)
    new_mat=np.subtract(mat1,matrix)
    return new_mat
         
def choose_neighbour(G,gamma,m,sim,t):
    dist=new_hyperbolic_dist(G,gamma,sim,t)
    connect_to = np.random.choice(G.nodes(),m, replace=False, p=dist)
    return(connect_to)


def new_hyperbolic_dist(G,gamma,sim,t):
    distance=[]
    beta = 1 / (gamma - 1)
    for node in G.nodes():
        rt = log(t)
        rs = (beta * log(int(node)) + (1 - beta) * (log(t)))
        theta =(sim[int(node) - 1][t - 1] * (pi))
        d1=cosh(2*rs)* cosh(2*rt)
        d2 = sinh(2 * rs) * sinh(2 * rt)*cos(theta)
        d=d1-d2
        if (d<1):  #due to precision problems, d<1 is set to 1 to get correct final h_distance, archcos=0
            d=1
            d=0.5*acosh(d)
            distance.append(d)
        else:
            d = acosh(d)
            distance.append(d)
    sum_dist = sum(distance)
    np.asarray(distance)
    np.nan_to_num(distance)# if nan or inf set to zero
    distance = [i / sum_dist for i in distance]
    return distance


#save degree in a text file for further analysis 
def gen_degree(G):
    plt.close() 
    degfile = open('degreelist.csv','w')
    deg_list=[]
    for n in G.nodes():
        deg_list.append(G.degree(n))
    degrees=np.array(deg_list)
    print(degrees)
    for d in deg_list:
        degfile.write(str(d)+'\n')
    return(degrees)


#powerlaw fitting and degree exponent.
def powerlaww(G):
    degree_list=gen_degree(G)
    fit=powerlaw.Fit(degree_list,linear_bins=True)
    print("The degree of exponent is:", fit.power_law.alpha)
    print("gamma's statndard error is:",fit.power_law.sigma)
    print("Likelihood of PL and Exponetial:",fit.distribution_compare('power_law','exponential'))
    fig1=fit.plot_pdf(color='b', linewidth=0, marker='o')
    #fig2=fit.plot_ccdf(color='b', linewidth=2)
    #fit.power_law.plot_pdf(color='b',linestyle='--')
    #fit.plot_cdf(color='r', linewidth=1)
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
    deg, =plt.plot(x,y,label='Degree Distribution',linewidth=0 , marker='o', markersize=8, color='r',alpha=alpha)
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

    plt.ylabel('P(k)')
    plt.xlabel('k')
    plt.title( 'Degree Distribution (Choice Method)')
    plt.show()


if __name__ == "__main__":
    Pop_Sim(4,3)


    
