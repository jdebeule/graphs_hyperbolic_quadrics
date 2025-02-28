#Initialisation
LoadPackage("fining");
q := 2;
n := 3;

#Setup of the graph G_n. We need Q+(2n+1,q), one generators and the points not in the generator.
ps := HyperbolicQuadric(2*n+1,q);
pg := AmbientSpace(ps);
generators := Solids(ps);
pi1 := First(generators,x->x=x); #Just take the first generator in the list.
pspts := AsList(Points(ps));;

#The vertices are the points of Q+(2n+1,q) not in pi1
vertices := Filtered(pspts,x->not x in pi1);;
group := CollineationGroup(ps);
stab := FiningStabiliser(group,pi1); #This will be a group stabilising vertices and respecting the adjacency.

#Adjacency function. Two different vertices are adjacent if they span a secant or (they are collinear and the connecting line meets pi_1)
adj := function(x,y)
if x<>y then
    if not IsCollinear(ps,x,y) then
        return true;
    else
        return ProjectiveDimension(Meet(Span(x,y),pi1))=0;
    fi;
else
    return false;
fi;
end;

#We construct the graph G_n with the usual grape command
g_n := Graph(stab,vertices,\^,adj,true);

#The graph NOplus(2n+2,q)
#We neede the same Q+(2n+1,q), but now the points of the ambient space not in the quadric.

pts := AsList(Points(pg));
verticesnoplus := Filtered(pts,x->not x in ps);

#The adjacency can be computed using the polarity associated to the quadric.
#Two vertices x,y are adjacent if the connecting line is a tangent line to
#the quadric, hence if y^perp contains x

delta := PolarityOfProjectiveSpace(ps);
adjnoplus := function(x,y)
if x = y then return false;
else return x in y^delta;
fi;
end;

#We construct the graph noplus with the usual grape command.
noplus := Graph(group,verticesnoplus,OnProjSubspaces,adjnoplus,true);

#The following can be slow for large graphs. This is just a nice to have and not necessary right now.
group1 := AutomorphismGroup(g_n); #can be slow for large graphs...
group2 := AutomorphismGroup(noplus);

#one way to show that g_n is not isomorphic with noplus is to investigate the cliques.

cliques_g_n := CompleteSubgraphs(g_n,-1,2);
cliques_noplus := CompleteSubgraphs(noplus,-1,2); #The result is clear: the graphs are non-isomorphic.

#It will be useful for the sequel to have adj1.
#This is a part of the adjacency relation of g_n, i.e. the vertices
#are adjacent if the connecting line is a secant.
adj1 := function(x,y)
if x<>y then
    if not IsCollinear(ps,x,y) then
        return true;
    else
        return false;
    fi;
else
    return false;
fi;
end;

#We have seen cliques of size 5 and 8

cliques5 := CompleteSubgraphsOfGivenSize(g_n,5,2,true);
cliques8 := CompleteSubgraphsOfGivenSize(g_n,8,2,true);

#We do the usual thing. Consider a clique as point set in the quadric and look at its span.
c5_1 := cliques5[1];
c5_1pts := VertexNames(g_n){c5_1};
space := Span(c5_1pts);
TypeOfSubspace(ps,space); #The result is clear, a maximal clisue of size 5 is an ellitpic quadric.

#look at the spans, but first note that a generator meeting pi1 in a plane, is indeed a clique of size 8. So
#we expect to find at least one generator in the spans.

c8_1 := cliques8[1];
c8_1pts := VertexNames(g_n){c8_1};
c8_2 := cliques8[2];
c8_2pts := VertexNames(g_n){c8_2};
Span(c8_2pts);
c8_3 := cliques8[3];
c8_3pts := VertexNames(g_n){c8_3};
Span(c8_3pts);
c8_4 := cliques8[4];
c8_4pts := VertexNames(g_n){c8_4};
Span(c8_4pts);

spaces8 := List(cliques8,x->Span(VertexNames(g_n){x}));

solids := Filtered(spaces8,x->ProjectiveDimension(x)=3);
List(solids,x->x in ps); #indeed, one solid is a generator!

t := Filtered(solids,x->not x in ps);
solid := t[1];
TypeOfSubspace(ps,solid);

#We have to find out how this solid intersects the quadric.
#We will use the adj1 function
#make sure you have the solid that is not a generator, seems to
#be c8_3pts.

Span(c8_2pts) in ps;
ProjectiveDimension(Span(c8_2pts));
adjacencies := List(c8_2pts,x->List(c8_2pts,y->adj1(x,y))); #We see that each point of the clique is non-adj1 with three others, hence, these adjacencies are of type 2, i.e. connecting lines
#are contained in Q+ and meet pi1. Hence the clique conists of the vertices in two planes meeting pi1 in the same line.

#this is the complete space, always interesting :-)
Span(c8_3pts)
adjacencies := List(c8_3pts,x->List(c8_3pts,y->adj1(x,y))); #We see only adjacencies1, hence these are 8 points of an ovoid!

#one more
space := Span(c8_1pts);
adjacencies := List(c8_1pts,x->List(c8_1pts,y->adj1(x,y)));
Meet(space,pi1); #each is adj1 with exactly 6 other points, and connects with exactly one point to a line meeting pi1.

#Let's compute all connecting lines and see where they meet.
pairs := List(c8_1pts,x->Filtered(c8_1pts,y->not adj1(x,y))); #This will include the point x itself, which is convenient, since then we can use the resulting pair to compute the connecting line
lines := List(pairs,x->Span(x));
List(lines,x->x in ps); #double check, should all be lines in ps.

Meet(lines); #seems to be a point!
point := Meet(lines);
point in pi1; #Now we know enough -> space is a cone with vertex "point" and base a Q-(3,q), meeting pi1 in a line. Indeed, we then see exactly 8 vertices of g_n.

#As an exercise, we can now see how many copies of each clique of size 5 and 8 we find, using the automorphim group of g_n.

o1 := Orbit(group1,Set(c8_1),OnSets);;
o2 := Orbit(group1,Set(c8_2),OnSets);;
o3 := Orbit(group1,Set(c8_3),OnSets);;
o4 := Orbit(group1,Set(c8_4),OnSets);;

Length(o1);
Length(o2);
Length(o3);
Length(o4);

o5 := Orbit(group1,Set(c5_1),OnSets);; #takes 10-20 seconds
Length(o5);

