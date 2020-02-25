We can find some info here: http://preprocessed-connectomes-project.org/abide/index.html

The dataset contains the data of 49 participants affected by Autism (a_1) Age [5,10)  
and 52 healthy controls (h_1) Age [5,10)  that underwent high-resolution structural magnetic resonance neuroimaging.
	
We build the brain network using the time series obtained from the fMRI. We compute pairwise correlation 
between these points and we obtain a complete graph $G=(V, E, p)$ with $|V|=116$. 
 

The weights (distances between nodes for us) for the graphs are compute using the coordinates of the nodes as 
given in the first three columns of the file Node_AAL116.node
therefore they are the same for all the graphs. The probabilities instead are obtained from the absolute
 value of the correlation coefficients.  

File name: |type-of-patiet a_1/h_1|_|ID-of-the-subject|.txt
