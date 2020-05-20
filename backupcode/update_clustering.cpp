#include <Rcpp.h>
using namespace Rcpp;

// Function to update the partions 
// given the old partition and a record we inspect, when its label switches 
// from oldlabel to newlabel, update the partition accordingly
//  clustering returns 
//  * clize: a list of n values representing cluster sizes
//  * ksize: a number of non-empty clusters
//  * clmembers: an nxn matrix, where each row corresponds to a cluster
//    and each row 'i' contains, from left to right, the indices which(lambda == i)
//    (with C convention, so indices start at zero)
//    ... the rest of the matrix 'clmembers' is filled with '-1'

// [[Rcpp::export]]


List update_clustering(const List clustering_old,
	const int j,
	const int oldlabel,
	const int newlabel){


	IntegerMatrix clmembers_new = clone(wrap(clustering_old["clmembers"]));
	IntegerVector ksize_vec = clone(wrap(clustering_old["ksize"]));
	IntegerVector clsize_new = clone(wrap(clustering_old["clsize"]));
	int ksize_new = ksize_vec(0);

	// Rprintf("looking at record label %i", j + 1);

	if (oldlabel == newlabel){
		// do nothing
	}else{
		// remove the jth record from the current partiton
		bool shift = false;
		for (int icol = 0; icol < clsize_new[oldlabel]; icol++){
			// move every record in the same cluster forward by one 
			if (clmembers_new(oldlabel, icol) == j){
				shift = true; // find the position of the jth record in the matrix clmembers
			}
			if (shift){
				// Rprintf("move record %i to the postion of record %i",clmembers_new(oldlabel, icol + 1), clmembers_new(oldlabel, icol) );
				// move next label to current place, and set the last one to -1
				if (icol + 1 < clsize_new[oldlabel]){
					clmembers_new(oldlabel, icol) = clmembers_new(oldlabel, icol + 1);
				}else{
					clmembers_new(oldlabel, icol) = - 1 ;
				}
			}
		}
		clsize_new[oldlabel]--;
		// if cluster becomes empty
		if (clsize_new[oldlabel] == 0){
			ksize_new --;
		}
		// add the jth record to the new cluster
		clsize_new[newlabel]++;
		clmembers_new(newlabel,clsize_new[newlabel] - 1) = j; 
		// if this creates. a new cluster
		if (clsize_new[newlabel] == 1){
			ksize_new++;
		}
	}

	return  List::create(Named("clsize") = clsize_new, Named("ksize") = ksize_new, Named("clmembers") = clmembers_new);
}