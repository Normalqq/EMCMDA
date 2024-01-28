# -*- codeing = utf-8 -*-
# @File : EMCMDA.py
# @Software: PyCharm
import math
import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc
import copy
import numpy.matlib





def SVT(Y, b, ty):
    #奇异值分解
    U,S,V= np.linalg.svd(Y)
    #奇异值收缩
    for index in range(0,S.size):#s.size函数是矩阵元素的个数
        if S[index] >= b*ty[index]:
            S[index] = S[index] - b*ty[index]
        else:
            S[index] = 0
    #奇异值矩阵标准化
    s = np.diag(S)
    row , col = Y.shape[0] , Y.shape[1]
    if row < col:
        s_n = np.column_stack((s, np.zeros((row, col - row))))
    else:
        s_n = np.row_stack((s, np.zeros((row-col, col))))
    #U*奇异值收缩矩阵*V
    Y_n = np.dot(U, np.dot(s_n, V))
    return Y_n






def MMTSPN(alpha, beta,gamma, t, omega, tol1, tol2, maxiter,A,B,p,r):
    X = t
    W = X
    Y = X
    iter0 = 1
    stop1 = 1
    stop2 = 1

    # the processing of computing Wt
    H=np.dot(B.transpose(),A)
    U,S,V= np.linalg.svd(H)#这里面奇异值的个数只有一个
    u,s,v= np.linalg.svd(X)
    W=np.zeros_like(s)
    for i in range(0,S.size):
        W[i]=p*(1-S[i])*pow(s[i],p-1)
    if(W[i]<0):
        W[i]=0

    while stop1 > tol1 or stop2 > tol2:
        # the processing of computing T
        tran = (1/beta) * (Y+alpha*(t*omega))+X
        T = tran - (alpha/(alpha+beta))*omega*tran
        T[T < 0] = 0
        T[T > 1] = 1

        # the processing of computing X
        X_1 = SVT(T-(1/beta)*E, 1/beta,W)

        # the processing of computing E
        E = E + gamma*(X_1-T)

        stop1_0 = stop1
        if np.linalg.norm(X) != 0:
            stop1 = np.linalg.norm(X_1-X) / np.linalg.norm(X)
        else:
            stop1 = np.linalg.norm(X_1-X)
        stop2 = np.abs(stop1-stop1_0)/(max(1, np.abs(stop1_0)))
        X = X_1

        if iter0 >= maxiter:
            iter0 = maxiter
            print('reach maximum iteration,did not converge!')
            break
        iter0 = iter0 + 1
    T_recover = T
    return T_recover, iter0


def run_MC(t):
    # MMTSPN的参数
    maxiter = 300
    alpha = 20
    beta = 5
    gamma = 1
    p=1
    tol1 = 2 * 1e-3
    tol2 = 1 * 1e-5
    omega = np.zeros(t.shape)
    omega[t.nonzero()] = 1
    #插入第一层循环，或者叫第一步
    for i in range(0,1):
        U, S, V = np.linalg.svd(t)
        r = 5
        A = U[:r, :]
        #print(np.dot(A,A.transpose()))
        B = V[:r, :]
        #print(np.dot(B,B.transpose()))

        t, k = MMTSPN(alpha, beta,gamma, t, omega, tol1, tol2, maxiter,A,B,p,r)
    Smmi = t
    H_1 = t_new[0:miRNA.shape[0], miRNA.shape[0]:T.shape[1]]
    return H_1

#T2.0为异构图关联矩阵，左上为miRNA相似性矩阵，右上为miRNA—疾病关联矩阵，右下为疾病相似性矩阵，左下为miRNA—疾病关联矩阵的转置,维数为878*878
#miRNA为miRNA相似性矩阵（495*495），疾病为疾病相似性矩阵（383*383），SM_miRNA_k为已知关联（5430*2），SM_miRNA_uk为未知关联（184155*2）
H = np.loadtxt(r'T2.0.txt', dtype=float)
miRNA = np.loadtxt(r"GIP_miRNA.txt",dtype=float)



if __name__ == "__main__":
    Scores_H = run_MC(H) #Scores_M是预测评分矩阵






