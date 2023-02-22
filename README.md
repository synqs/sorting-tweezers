# sorting-tweezers
Code to assemble a defect free array of tweezers.

For FIFO mode, the method described here allows to use very high sample rates from the AWG even when the computer cannot calculate data very fast. The overall scheme involves using two data arrays in computer memory which are used alternately to stream their contents to the AWG. The selection of which array will stream to the AWG is done using a Boolean variable. When one array is being used for streaming to the AWG, the other array is used for calculating and storing the data of the next move. After this calculation, the Boolean variable is flipped, and the role of the two arrays interchange (i.e. the first one is now used for calculation and the second one for streaming). Since the array used for streaming to the AWG is not used for calculating, this scheme guarantees that there will never be a shortage of data when the AWG is outputting at very high sample rates.

For more details see PhD thesis of Rohit Prasad Bhatt (publicly available at https://doi.org/10.11588/heidok.00031968).
