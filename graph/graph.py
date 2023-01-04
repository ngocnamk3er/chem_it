import matplotlib.pyplot as plt
import numpy as np

class stage:
    def __init__(self,label,w):
        self.label = label
        self.w = w
        self.nextstage = list()

ra = stage("ra",0.0)
comp1 = stage("comp1",-13.1)
comp2 = stage("comp2",-3.7)
comp3 = stage("comp3",-12.8)
ts1 = stage("ts1",-4.7)
ts2 = stage("ts2",-6.0)
ts3 = stage("ts3",20.4)
ts4 = stage("ts4",0.0)
ts5 = stage("ts5",0.0)
ts6 = stage("ts6",0.0)
ts7 = stage("ts7",0.0)
ts8 = stage("ts8",0.0)
ts9 = stage("ts9",0.0)
ts10 = stage("ts10",0.0)
ts11 = stage("ts11",0.0)
ts12 = stage("ts12",0.0)
ts13 = stage("ts13",0.0)
is1  = stage("is1",0.0)
is2 = stage("is2",0.0)
is3 = stage("is3",0.0)
is4 = stage("is4",0.0)
pr1 = stage("pr1",0.0)
pr2 = stage("pr2",0.0)
pr3 = stage("pr3",0.0)
pr4 = stage("pr4",0.0)
pr5 = stage("pr5",0.0)

# make tree

ra.nextstage.append(comp1)
ra.nextstage.append(comp2)
ra.nextstage.append(comp3)
comp1.nextstage.append(ts1)
comp1.nextstage.append(ts2)
comp1.nextstage.append(ts3)



def show_graph(root_stage):

    wArr = list()
    labelArr = list()
    def tree_bfs(root_stage):
        visited = list()
        queue = list()
        queue.append(root_stage)
        visited.append(root_stage.label) 
        wArr.append(root_stage.w)
        labelArr.append(root_stage.label)
        while queue:
            s = queue.pop()
            print(s.label)
            for i in s.nextstage:
                queue.append(i)
                wArr.append(i.w)
                labelArr.append(i.label)
    tree_bfs(root_stage)
    label1 = np.array(labelArr)
    w1 = np.array(wArr)
    plt.plot( label1,w1)
    plt.plot( np.array(labelArr),np.array([1,2,3,4,5,6,7]))
    plt.show()

def main(): 
    show_graph(ra)

if __name__ == "__main__":
    main()