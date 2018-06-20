The conf directories here are named with a few things in mind.

-Number of years of data. In the event of 1yr, the conf will simply state the name of that year.
-Splined IC79. For such analyses (later ones, my tests and analysis) I will include the letter 'b' at the end of the conf. If I use sirin's IC79, nothing.
-Whose data: if I'm using stefan's data (usually what I'm doing), I won't add anything. If I'm using data that I myself pull corrected, I'll put a 'm' at the end. The script will know to do with because of a conf variable, mydata.
-What binning: if I'm using stefan's binning (usually what I'm doing), I'll use a conf variable, binning, which can be npz_restrict, npz_full, uniform. Uniform conf has a 'u' in it.
-Finally, I want to know when I look at the conf folder what it is being used for. I will include a README.txt file in each conf to detail this. 
-One more thing - 7yr is assumed to contain MESE data. if I want there _not_ to be MESE data, I'll put 'noMESE' at the end.

The most important thing is that I know EXACTLY which conf was used for which thing. If I use one for something, I have to make sure to write it in the appropriate README file. If I no longer use a particular conf for something, I have to delete it from the README AND from the data storage.
