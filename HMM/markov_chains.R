#!/usr/bin/env Rscript
#markov_chains.R
#!/usr/bin/env Rscript
# =============================================================================
#
#         FILE: markov_chains.R
#
#        USAGE: Rscript markov_chains.R
#
#  DESCRIPTION: This program show you how to generate a DNA sequence using a 
#               multinomial model of DNA sequence evolution, and then moves onto
#               a Markov model of DNA sequence evolution, plotting data as we 
#               go along
#
#       AUTHOR: This program has parts from Benjamin Tovar and Avril Coghlan, 
#               and Chesley Leslin
#
# =============================================================================

# We wil use the "seqinr" package and other packages
# For the installation of the package, type:
# install.packages("seqinr")
library(plyr)
library(seqinr)
library(ggplot2)
library(gridExtra)

# Create a function to plot the first order markov data
plotData <- function(seq,
                                          nucleotides,
                                                               title){
      #    seq              = sequence generated
      #    nucleotides      = alphabet that compounds the sequence 
      #    tittle           = title of the plot

      # ***** Study the composition bias of the sequence *****
      # Count the frequency of each base 
      # in the sequence using the "count" function
      freqNucleotide <- count(seq,1,alphabet=nucleotides,freq=TRUE)

  # Count the frequency of dinucleotides 
  # in the sequence using the "count" function
  freqDiNucleotide <- count(seq,2,alphabet=nucleotides,freq=TRUE)

    # Count the frequency of trinucleotides 
    # in the sequence using the "count" function
    freqTriNucleotide <- count(seq,3,alphabet=nucleotides,freq=TRUE)

    # We have to use data frames in ggplot 2, so just covert, and the give column names
    df <- as.data.frame(freqNucleotide)
      colnames(df) <- c("Base", "Base_Proportion")
      p1 <- ggplot(df, aes(x = Base, y = Base_Proportion, fill=Base )) +
          geom_bar(stat = "identity")
        p1 <- p1 + theme(legend.position="none") +
            ggtitle("Compositional bias of each nucleotide")

          df <- as.data.frame(freqDiNucleotide)
            colnames(df) <- c("Base", "Base_Proportion")
            p2 <- ggplot(df, aes(x = Base, y = Base_Proportion, fill=Base )) +
                geom_bar(stat = "identity")
              p2 <- p2 + theme(legend.position="none") +
                  ggtitle("Compositional bias of each dinucleotide")

                df <- as.data.frame(freqTriNucleotide)
                  colnames(df) <- c("Base", "Base_Proportion")
                  p3 <- ggplot(df, aes(x = Base, y = Base_Proportion, fill=Base )) +
                      geom_bar(stat = "identity")
                    p3 <- p3 + theme(legend.position="none") +
                        ggtitle("Compositional bias of each trinucleotide") +
                            theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=8) )

                      # Now, plot the results in the same plot:
                      grid.arrange(p1, p2, p3, nrow=3, top=title)
}

# Create a function to generate the sequence.
generateFirstOrderSeq <- function(lengthSeq,
                                                                nucleotides,
                                                                                              initialProb,
                                                                                                                            firstOrderMatrix){
      #    lengthSeq        = length of the sequence
      #    nucleotides      = alphabet that compounds the sequence 
      #    initialProb      = initial probability distribution
      #    firstOrderMatrix = matrix that stores the probability distribution of a 
      #                       first order Markov Chain

      # Create a vector for storing the new sequence
      outputSeq <- character()
  firstnucleotide <- sample(nucleotides, 1, rep=TRUE, prob=initialProb)
    # Store the nucleotide for the first position of the sequence
    outputSeq[1]    <- firstnucleotide

    # Let the computer decide:
    for(i in 2:lengthSeq){
            prevNuc <- outputSeq[i-1]
        currentProb <- firstOrderMatrix[prevNuc,]
            #cat("\n", prevNuc, " " , currentProb, sep=" ")
            outputSeq[i] <- sample(nucleotides,1,prob=currentProb)
          }
      return(outputSeq)
}
# Function to generate a DNA sequence, given a HMM and the length of the sequence to be generated.
generateFirstOrderhmmseq <- function(lengthSeq, nucleotides, initialProb, states, transitionmatrix, emissionmatrix){

      #    lengthSeq        = length of the sequence
      #    nucleotides      = alphabet that compounds the sequence 
      #    initialProb      = initial probability distribution
      #    states           = states for the HMM
      #    transitionmatrix = matrix that stores the probability distribution of transitions
      #    emissionmatrix   = matrix that stores the probability distribution of emissions

      outputSeq      <- character()             # Create a vector for storing the new sequence
  mystates        <- character()             # Create a vector for storing the state that each position in the new sequence
    # was generated by
    # Choose the state for the first position in the sequence:
    firststate      <- sample(states, 1, rep=TRUE, prob=initialProb)
    # Get the probabilities of the current nucleotide, given that we are in the state "firststate":
    probabilities   <- emissionmatrix[firststate,]
      # Choose the nucleotide for the first position in the sequence:
      firstnucleotide <- sample(nucleotides, 1, rep=TRUE, prob=probabilities)
      outputSeq[1]   <- firstnucleotide         # Store the nucleotide for the first position of the sequence
        mystates[1]     <- firststate             # Store the state that the first position in the sequence was generated by

        for (i in 2:lengthSeq){
                prevstate    <- mystates[i-1]           # Get the state that the previous nucleotide in the sequence was generated by
            # Get the probabilities of the current state, given that the previous nucleotide was generated by state "prevstate"
            stateprobs   <- transitionmatrix[prevstate,]
                # Choose the state for the ith position in the sequence:
                state        <- sample(states, 1, rep=TRUE, prob=stateprobs)
                # Get the probabilities of the current nucleotide, given that we are in the state "state":
                probabilities <- emissionmatrix[state,]
                    # Choose the nucleotide for the ith position in the sequence:
                    nucleotide   <- sample(nucleotides, 1, rep=TRUE, prob=probabilities)
                    outputSeq[i] <- nucleotide             # Store the nucleotide for the current position of the sequence
                        mystates[i]  <- state                  # Store the state that the current position in the sequence was generated by
                      }

          for (i in 1:lengthSeq){
                  nucleotide   <- outputSeq[i]
              state        <- mystates[i]
                  print(paste("Position", i, ", State", state, ", Nucleotide = ", nucleotide))
                }
          return(outputSeq)
}

viterbi <- function(sequence, transitionmatrix, emissionmatrix){

      #    sequence         = character vector for the sequence
      #    transitionmatrix = matrix that stores the probability distribution of transitions
      #    emissionmatrix   = matrix that stores the probability distribution of emissions

      # This carries out the Viterbi algorithm.
      # Adapted from "Applied Statistics for Bioinformatics using R" by Wim P. Krijnen, page 209
      # ( cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf )
      # Get the names of the states in the HMM:
      states <- rownames(emissionmatrix)

  # Make the Viterbi matrix v:
  v <- makeViterbimat(sequence, transitionmatrix, emissionmatrix)

    # Go through each of the rows of the matrix v (where each row represents
    # a position in the DNA sequence), and find out which column has the
    # maximum value for that row (where each column represents one state of
    # the HMM):
    mostprobablestatepath <- apply(v, 1, function(x) which.max(x))

    # Print out the most probable state path:
    prevnucleotide            <- sequence[1]
      prevmostprobablestate     <- mostprobablestatepath[1]
      prevmostprobablestatename <- states[prevmostprobablestate]
        startpos <- 1
        for (i in 2:length(sequence)){
                nucleotide            <- sequence[i]
            mostprobablestate     <- mostprobablestatepath[i]
                mostprobablestatename <- states[mostprobablestate]
                if (mostprobablestatename != prevmostprobablestatename){
                          print(paste("Positions",startpos,"-",(i-1),
                                                        "Most probable state = ", prevmostprobablestatename))
                      startpos <- i
                          }
                    prevnucleotide <- nucleotide
                        prevmostprobablestatename <- mostprobablestatename
                      }
          print(paste("Positions",startpos,"-",i,
                                    "Most probable state = ", prevmostprobablestatename))
}

makeViterbimat <- function(sequence, transitionmatrix, emissionmatrix) {

      #    sequence         = character vector for the sequence
      #    transitionmatrix = matrix that stores the probability distribution of transitions
      #    emissionmatrix   = matrix that stores the probability distribution of emissions  

      # This makes the matrix v using the Viterbi algorithm.
      # Adapted from "Applied Statistics for Bioinformatics using R" by Wim P. Krijnen, page 209
      # ( cran.r-project.org/doc/contrib/Krijnen-IntroBioInfStatistics.pdf )
      # Change the sequence to uppercase
      sequence <- toupper(sequence)
  # Find out how many states are in the HMM
  numstates <- dim(transitionmatrix)[1]
    # Make a matrix with as many rows as positions in the sequence, and as many
    # columns as states in the HMM
    v <- matrix(NA, nrow = length(sequence), ncol = dim(transitionmatrix)[1])
    # Set the values in the first row of matrix v 
    # (representing the first position of the sequence) to 0
    v[1, ] <- 0
      # Set the value in the first row of matrix v, first column to 1
      v[1,1] <- 1
      # Fill in the matrix v:
      for (i in 2:length(sequence)) { # For each position in the DNA sequence:
              for (l in 1:numstates) { # For each of the states of in the HMM:
                        # Find the probabilility, if we are in state l, 
                        # of choosing the nucleotide at position in the sequence
                        statelprobnucleotidei <- emissionmatrix[l,sequence[i]]

            # v[(i-1),] gives the values of v for the (i-1)th row of v, 
            # ie. the (i-1)th position in the sequence.
            # In v[(i-1),] there are values of v at the (i-1)th 
            # row of the sequence for each possible state k.
            # v[(i-1),k] gives the value of v at the (i-1)th row 
            # of the sequence for a particular state k.

            # transitionmatrix[l,] gives the values in the lth row of the 
            # transition matrix, xx should not be transitionmatrix[,l]?
            # probabilities of changing from a previous state k to a current state l.

            # max(v[(i-1),] * transitionmatrix[l,]) 
            # is the maximum probability for the nucleotide observed
            # at the previous position in the sequence in state k, 
            # followed by a transition from previous
            # state k to current state l at the current nucleotide position.

            # Set the value in matrix v for row i (nucleotide position i), 
            # column l (state l) to be:
            v[i,l] <-  statelprobnucleotidei * max(v[(i-1),] * transitionmatrix[,l])
                }
        }
        return(v)
}

#############End of Functions

pdf("markov_plots.pdf")

#    Example of a multinomial model of DNA evoloution (the current nucleotide is
#    totally "independent" of the previous nucleotide).

#    The multinomial model is:
#    p(A)+p(C)+p(G)+p(T) = 1.0
#    0.2 +0.3 +0.3 +0.2  = 1.0

# Define the DNA alphabet that will be used to put names to the objects
nucleotides <- c("A","C","G","T")
# Create the vector that represents the probability distribution of the model
zeroOrderProbablities <- c(0.2,0.3,0.3,0.2)
# Put the name of reference of each base 
names(zeroOrderProbablities) <- nucleotides
# Create a sequence of 1000 bases using this model.
zeroOrderSeq <- sample(nucleotides,1000,rep=T,prob=zeroOrderProbablities)

plotData(zeroOrderSeq, nucleotides, "Multinomial Model of DNA Evoloution")


#    Example of a Markov Chain of first order (the current nucleotide only
#    depends on the previous nucleotide).

#    The multinomial model per base is:
#    p(A|A)+p(C|A)+p(G|A)+p(T|A) = 1.0
#    p(A|C)+p(C|C)+p(G|C)+p(T|C) = 1.0
#    p(A|G)+p(C|G)+p(G|G)+p(T|G) = 1.0
#    p(A|T)+p(C|T)+p(G|T)+p(T|T) = 1.0        

#    So:
#    0.2  + 0.3   + 0.3   + 0.2   = 1.0
#    0.1  + 0.41  + 0.39  + 0.1   = 1.0
#    0.25 + 0.25  + 0.25  + 0.25  = 1.0
#    0.5  + 0.17  + 0.17  + 0.17  = 1.0

# Add the probability distribution per base:
afterAprobs <- c(0.2,0.3,0.3,0.2)  # Set the values of the probabilities, where the previous nucleotide was "A"
afterCprobs <- c(0.1,0.41,0.39,0.1)  # Set the values of the probabilities, where the previous nucleotide was "C"
afterGprobs <- c(0.25,0.25,0.25,0.25)# Set the values of the probabilities, where the previous nucleotide was "G"
afterTprobs <- c(0.5,0.17,0.17,0.17)# Set the values of the probabilities, where the previous nucleotide was "T"

# Create the matrix that will store the probability distribution given 
# a certain nucleotide:
mytransitionmatrix <- matrix(c(afterAprobs, afterCprobs, afterGprobs, afterTprobs), 4, 4, byrow = TRUE) # Create a 4 x 4 matrix

# Put names to the 2 dimensions of the matrix 
colnames(mytransitionmatrix) <- nucleotides
rownames(mytransitionmatrix) <- nucleotides

# Now we got a matrix
#    > mytransitionmatrix
#         A    C    G    T
#    A 0.60 0.10 0.10 0.20
#    C 0.10 0.50 0.30 0.10
#    G 0.05 0.20 0.70 0.05
#    T 0.40 0.05 0.05 0.50

# In order to continue, we need an initial probability distribution to know
# which base is the most probable to start up the sequence.
inProb <- c(0.4,0.1,0.1,0.4)
names(inProb) <- nucleotides
# So, the sequence will have a 40% to start with an A or a T and 10% with C or G

# Use the generateFirstOrderSeq function to generate a sequence of 1000 bases long
firstOrderSeq <- generateFirstOrderSeq(1000,nucleotides,inProb,mytransitionmatrix)
# plot the data
plotData(firstOrderSeq, nucleotides, "Markov Chain of first order")

# Now let's create a much larger character vector, that we'll append to in a loop
#all <- character() # create an empty vector
#for(i in 1:100){
#  firstOrderSeq <- generateFirstOrderSeq(1000,nucleotides,inProb,mytransitionmatrix)
#  all <- c(firstOrderSeq, all)
#  #cat("Dna " , i , " " , firstOrderSeq, "\n", sep="")
#} 
#plotData(all, nucleotides,  "Markov Chain of first order with Expanded Sequence")

#    Example of a Hidden Markov Chain of first order (the current nucleotide only
#    depends on the previous nucleotide).

states              <- c("AT-rich", "GC-rich") # Define the names of the states
ATrichprobs         <- c(0.7, 0.3)             # Set the probabilities of switching states, where the previous state was "AT-rich"
GCrichprobs         <- c(0.1, 0.9)             # Set the probabilities of switching states, where the previous state was "GC-rich"
theTransitionMatrix <- matrix(c(ATrichprobs, GCrichprobs), 2, 2, byrow = TRUE) # Create a 2 x 2 matrix
rownames(theTransitionMatrix) <- states
colnames(theTransitionMatrix) <- states

ATrichstateprobs    <- c(0.39, 0.1, 0.1, 0.41) # Set the values of the probabilities, for the AT-rich state
GCrichstateprobs    <- c(0.1, 0.41, 0.39, 0.1) # Set the values of the probabilities, for the GC-rich state
theEmissionMatrix   <- matrix(c(ATrichstateprobs, GCrichstateprobs), 2, 4, byrow = TRUE) # Create a 2 x 4 matrix
rownames(theEmissionMatrix) <- states
colnames(theEmissionMatrix) <- nucleotides

# We can use the generateFirstOrderhmmseq() function to generate a sequence using a particular HMM. 
# For example, to create a sequence of 30 nucleotides using the HMM with “AT-rich” and “GC-rich” states described in the 
# transition matrix thetransitionmatrix, the emission matrix theemissionmatrix, 
# and uniform starting probabilities (ie. ΠAT-rich = 0.5, ΠGC-rich = 0.5), we type:
initialProb <- c(0.5, 0.5)
hmmfirstOrderSeq = generateFirstOrderhmmseq(10000, nucleotides, initialProb, states, theTransitionMatrix, theEmissionMatrix)
plotData(hmmfirstOrderSeq, nucleotides,  "Hidden Markov Model of first order")

# Given a HMM, and a particular DNA sequence, you can use the Viterbi function 
# to find the state of that HMM that was most likely to have generated the nucleotide at each 
# position in the DNA sequence:

myseq <- c("A", "A", "G", "C", "G", "T", "G", "G", "G", "G", "C", "C", "C", "C",
                      "G", "G", "C", "G", "A", "C", "A", "T", "G", "G", "G", "G", "T", "G",
                                 "T", "C")
viterbi(myseq, theTransitionMatrix, theEmissionMatrix)


# for more advanced Markov chain analysis see the package: 
# Package 'markovchain' - R Project

dev.off()
