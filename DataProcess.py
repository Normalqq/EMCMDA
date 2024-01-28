# -*- codeing = utf-8 -*-
# @File : DataProcess.py
# @Software: PyCharm
import math
import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc
import copy
import numpy.matlib

#加载数据资源
miRNA = np.loadtxt(r"miRNA相似性矩阵.txt",dtype=float)
disease = np.loadtxt(r"疾病相似性矩阵.txt",dtype=float)
miRNA_disease_M = np.loadtxt(r"miRNA-疾病关联矩阵.txt", dtype=int)
miRNA_disease_k = np.loadtxt(r"miRNA-疾病-已知关联.txt",dtype=int)
miRNA_disease_uk = np.loadtxt(r"miRNA-疾病-未知关联.txt",dtype=int)




#计算miRNA高斯轮廓核相似性
def Gaussian():
    row=495
    sum=0
    miRNA1=np.matlib.zeros((row,row))
    for i in range(0,row):
        a=np.linalg.norm(miRNA_disease_M[i,])*np.linalg.norm(miRNA_disease_M[i,])
        sum=sum+a
    ps=row/sum
    for i in range(0,row):
        for j in range(0,row):
            miRNA1[i,j]=math.exp(-ps*np.linalg.norm(miRNA_disease_M[i,]-miRNA_disease_M[j,])*np.linalg.norm(miRNA_disease_M[i,]-miRNA_disease_M[j,]))
            if(miRNA[i,j]==0):
                miRNA[i,j]=miRNA1[i,j]
    df = pd.DataFrame(miRNA)
    df.to_csv("./GIP_miRNA.csv")
    return miRNA
#计算disease高斯轮廓核相似性
def Gaussian1():
    column=383
    sum=0
    disease1=np.matlib.zeros((column,column))
    for i in range(0,column):
        a=np.linalg.norm(miRNA_disease_M[:,i])*np.linalg.norm(miRNA_disease_M[:,i])
        sum=sum+a
    ps=column/sum
    for i in range(0,column):
        for j in range(0,column):
            disease1[i,j]=math.exp(-ps*np.linalg.norm(miRNA_disease_M[:,i]-miRNA_disease_M[:,j])*np.linalg.norm(miRNA_disease_M[:,i]-miRNA_disease_M[:,j]))
            if(disease[i,j]==0):
                disease[i,j]=disease1[i,j]
    df = pd.DataFrame(disease)
    df.to_csv("./GIP_disease.csv")
    return disease

if __name__ == "__main__":
    Gaussian()
    Gaussian1()





