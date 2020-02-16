The dataset contains the data of 54 participants with bipolar disorder (BP) (ids from 1 to 55 but subject 17 is missing) 
and 65 healthy controls (HC) that underwent high-resolution structural magnetic resonance neuroimaging.
	
We build the brain network using the time series with 118 time points obtained from the fMRI. We compute pairwise correlation 
between these points and we obtain a complete graph $G=(V, E, p)$ with $|V|=116$. The set of subjects with Bipolar disorder can be divided 
in 2 subsets: 36 over 54 are lithium-treated patients (BP-L) while 17 were not treated with lithium (BP-nL), 
the data for subject 34 are missing.
 

The weights (distances between nodes for us) for the graphs are compute using the coordinates of the nodes as 
given in the first three columns of the file Node_AAL116.node
therefore they are the same for all the graphs. The probabilities instead are obtained from the absolute
 value of the correlation coefficients.  

File name: |type-of-patiet BP, BP-nL, HC|_sbj_|ID-of-the-subject|.txt
