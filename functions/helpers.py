#########################################
# Author: Chenfu Shi
# Email: chenfu.shi@postgrad.manchester.ac.uk


# helper functions that can be used in more advanced settings eg. specific plotting function or data management that can be used from any file
# not to be used for actual pipeline run  

#########################################


# examples:
# file opener function



# plotter of weird things such as
# matplotlib.pyplot.plot(numpy.array(sorted(CSR_mat.sum(axis=0).tolist()[0])).cumsum()/CSR_mat.sum())

# matplotlib.pyplot.plot(numpy.array(sorted(CSR_mat.diagonal().tolist())).cumsum()/CSR_mat.diagonal().sum())

# fig, ax = matplotlib.pyplot.subplots()

# lengths=[]
# for i in frag_prop:
    # lengths.append(i[3])
    # ax.hist(lengths, bins= 100,range=(0, 10000))
# fig.show()