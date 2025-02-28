# graphs_hyperbolic_quadrics
Code to deal with some graphs from hyperbolic quadrics

This repository contains some GAP code to deal with graphs coming from hyperbolic quadrics. This code has been used to investigate cliques of a particular graph in order to show that it is non isomorphic to a co-spectral graph. 

There are two graphs involved, (i) the graph G_n and (2), the graph NO+(2n+1,q), both for n=3. To described both graphs, first choose a hyperbolic quadric Q+(2n+1,q) in the ambient projective space PG(2n+1,q). 

(i) The graph G_n 

Choose a generator pi_1 of Q+(2n+1,q). The vertices are the points of Q+(2n+1,q) not in pi1. Two vertices are adjacent if (1) the connecting line is a secant to the quadric or (2) the connecting line is contained in the quadric and meets pi1 in a point. 

(ii) The graph NO+(2n+1,q). The vertices are the points of PG(2n+1,q) not in Q+(2n+1,q), two vertices are adjacent if the connecting line is tangent to the quadric. 

To show that both graphs are non-isomorhpic for n=3, we show that the maximal cliques are not the same. We also classify the maximal cliques of G_n by simply computing the orbit representatives. As an easy exercise, we also characerize the different maximal cliques geometrically. 
