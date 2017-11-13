#TODO: some visualization
#for (i in 1:6) {plot(medoids_ascii[i,1:200],type="l",ylim=c(-5,5),col=i);par(new=TRUE)}
#...
#PLOT:
#plot manifold 2D distances WER / --> récupérer les distances ? quand ?
#fenetre tempo forme des courbes / --> OK (jour[type] / semaine, indices en arg)
#medoids / --> OK (moyennés sur 1 jour / type de jour / semaine)
#gain en prevision: clust puis full --> enercast (comment l'utiliser ?)
#
#> plot(cr$medoids[1:100,1],type="l")
#> #for (i in 1:15)plot(cr$medoids[1:100,1],type="l")
#> r = range(cr$medoids[1:96,])
#> for (i in 1:15) {plot(cr$medoids[1:96,i],type="l",ylim=r,col=i); par(new=TRUE) }
#> 
#
