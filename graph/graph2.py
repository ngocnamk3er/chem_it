import matplotlib.pyplot as plt
import numpy as np

class stage:
    def __init__(self,label,w,time):
        self.label = label
        self.w = w
        self.time = time
        self.nextstage = list()

ra = stage("ra",0.0,1)
comp1 = stage("comp1",-13.1,2)
comp2 = stage("comp2",-3.7,3)
comp3 = stage("comp3",-12.8,4)
ts1 = stage("ts1",-4.7,5)
ts2 = stage("ts2",-6.0,6)
ts3 = stage("ts3",20.4,8)
ts4 = stage("ts4",-49.5,9)
ts5 = stage("ts5",-41.3,11)
ts6 = stage("ts6",-76.7,18)
ts7 = stage("ts7",-73.2,16)
ts8 = stage("ts8",-57.5,12)
ts9 = stage("ts9",-70.2,13)
ts10 = stage("ts10",-63.2,19)
ts11 = stage("ts11",-57.5,15)
ts12 = stage("ts12",-66.1,20)
ts13 = stage("ts13",-68.7,21)
is1  = stage("is1",-91.0,7)
is2 = stage("is2",-96.5,10)
is3 = stage("is3",-103.5,17)
is4 = stage("is4",-94.5,11)
pr1 = stage("pr1",13.1,22)
pr2 = stage("pr2",-84.3,23)
pr3 = stage("pr3",-80.9,25)
pr4 = stage("pr4",-70.2,24)
pr5 = stage("pr5",-68.5,26)

# make tree

ra.nextstage.append(comp1)
ra.nextstage.append(comp2)
ra.nextstage.append(comp3)

comp1.nextstage.append(ts1)
comp2.nextstage.append(ts2)
comp3.nextstage.append(ts3)

ts1.nextstage.append(is1)
ts2.nextstage.append(is2)
ts3.nextstage.append(is3)
ts4.nextstage.append(pr1)
ts4.nextstage.append(pr2)
ts5.nextstage.append(is4)
ts6.nextstage.append(pr2)
ts7.nextstage.append(pr3)
ts8.nextstage.append(is4)
ts9.nextstage.append(pr4)
ts10.nextstage.append(pr5)
ts11.nextstage.append(is3)
ts12.nextstage.append(pr5)
ts13.nextstage.append(pr3)

is1.nextstage.append(ts5)
is1.nextstage.append(ts4)

is2.nextstage.append(ts8)
is2.nextstage.append(ts9)
is2.nextstage.append(ts10)
is2.nextstage.append(ts11)

is3.nextstage.append(ts12)
is3.nextstage.append(ts13)

is4.nextstage.append(ts6)
is4.nextstage.append(ts7)


visited = dict()

listLabel = list()
listW = list()


def dfs_tree(root_stage,listLabel,listW):
    listLabel.append(root_stage.label)
    listW.append(root_stage.w)
    if (len(root_stage.nextstage)==0):
        # plt.xticks(np.array(listTime),np.array(listLabel))
        plt.plot(np.array(listLabel),np.array(listW))
        # plt.show()
    for i in root_stage.nextstage:
        # if i.label not in visited:
        #     visited[i.label] = 1
            dfs_tree(i,listLabel,listW)
            listLabel.pop()
            listW.pop()
        # else:
        #     print(i.label, "is already visited")

def show_graph(root_stage):
    
    pass


def main(): 
    dfs_tree(ra,listLabel,listW)
    plt.show()

if __name__ == "__main__":
    main()
    # x = np.array([4,1,2,3])
    # y = np.array([20,21,22,23])
    # my_xticks = ['John','Arnold','Mavis','Matt']
    # plt.xticks(x, my_xticks)
    # plt.plot(x, y)
    # plt.show()