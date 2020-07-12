import numpy as np
import pandas as pd
import random
import time
import matplotlib.pyplot as plt
"""
选取30天的数据
前10天按照算法规则初始化堆场
后20天按照优化后的算法规则进场/不优化的算法规则进场
"""
# 堆场内区段最大容量 17*15*3 = 765
BN = 17
LN = 15
TN = 3
MAX_CONTAIN = 765
# 堆场区段数量
DN = 20
# 分组数（目的港有10个. 泊位有7个）
GN = 70
# 阈值
alpha = 0.95
# 堆场内最大饱和度
theta = 0.85

data = pd.read_excel(r'.\data.xls')
data = data[['月份','合计箱量','内贸箱量','外贸箱量','20尺箱量', '40尺箱量']]
data.columns = ['month', 'total', 'import', 'export', '20num', '40num']
data['month'] = pd.to_datetime(data.month).dt.month
data['rate'] = data['import'] / data['total']
data['import20num'] = (data['20num'] * data.rate).map(int)
data['import40num'] = (data['40num'] * data.rate).map(int)
data = data.dropna(axis=0)

# 取1月数据
num = data['import'][0]
fakeData = pd.DataFrame({'targetPort':[], 'importDay':[], 'exportDay':[], 'berth':[]})
fakeData['targetPort'] = np.random.randint(1, 11, size=num) # 10个目的港
fakeData['importDay'] = np.random.randint(1, 31, size=num) # 交箱时间
fakeData['exportDay'] = np.random.randint(1, 31, size=num) # 装船时间
fakeData['berth'] = np.random.randint(1, 8, size=num) # 7个泊位
fakeData['weight'] = np.random.randint(20, 28, size=num) # 随机重量级
fakeData['group'] = [0 for _ in range(num)]
fakeData.to_csv('fakeData.csv', index=None)
# fakeData = pd.read_csv('fakeData.csv')
# print(len(fakeData[fakeData.importDay > fakeData.exportDay]))
# print(len(fakeData[fakeData.importDay <= fakeData.exportDay]))
print(fakeData)
# # TODO 设置目的港偏好

def initData(fakeData):
    group = dict()
    for index in fakeData.index:
        item = (fakeData.targetPort[index], fakeData.importDay[index], fakeData.exportDay[index], fakeData.berth[index])
        ii = item[0]
        jj = item[3]
        # print('targetPort:', ii, 'berth:', jj)
        # targetPort
        for i in range(1, 11):
            if ii == i:
                # berth
                for j in range(1, 8):
                    if jj == j:
                        if ii * 10 + jj not in group:
                            group[ii * 10 + jj] = list()
                        group[ii * 10 + jj].append(index)
                        # print(index, 'assign', ii * 10 + jj)
                        break
                break
    # 由组内箱量降序排序，GN=70
    group = sorted(group.items(), key=lambda x:len(x[1]), reverse=True)
    for item in group:
        for index in item[1]:
            fakeData['group'][index] = item[0]
    fakeData.to_csv('fakeData.csv', index=None)

def initDock(boxStream):
    initDock = dict()
    for n in range(DN):
        dock = list()
        for b in range(BN):
            dock.append(list())
            for _ in range(LN):
                dock[b].append(list())
        initDock[n + 1] = dock
    # initd, initr = myRules(initDock, boxStream)
    # print('distance:', initd, 'reshuffle:',initr)
    return initDock
    
def distance(initDock):
    # distance
    dtmp = list()
    countBox = 0
    for nn in range(DN):
        dock = initDock[nn + 1]
        for b in range(BN):
            for l in range(LN):
                boxTie = dock[b][l]
                for box in boxTie:
                    countBox = countBox + 1
                    # countBoxt = 0
                    dgtmp = list()
                    for nnt in range(DN):
                        dockt = initDock[nnt + 1]
                        for bt in range(BN):
                            for lt in range(LN):
                                # 同场区，同贝位，同排位，一共1-3个箱子, 加上也无妨。
                                if nnt == nn and bt == b and lt == l:
                                    continue
                                boxTiet = dockt[bt][lt]
                                for boxt in boxTiet:
                                    # countBoxt = countBoxt + 1
                                    # 同组
                                    if box[3] == boxt[3]:
                                        dgtmp.append(abs(box[1] - boxt[1]))
                    dtmp.append(sum(dgtmp))
                # print(countBoxt)
    d = sum(dtmp)
    print(d)
    return d

def reshuffle(initDock):
    # reshuffle
    rtmp = list()
    countBox = 0
    for nn in range(DN):
        dock = initDock[nn + 1]
        for b in range(BN):
            for l in range(LN):
                rbltmp = 0
                boxTie = dock[b][l]
                bexportDays = list()
                bindexs = list()
                for index, box in enumerate(boxTie):
                    bexportDays.append(box[6])
                    bindexs.append(index)
                    for bexportDay, bindex in zip(bexportDays, bindexs):
                        if box[6] > bexportDay:
                            rbltmp = rbltmp + index - bindex
                            break
                    countBox = countBox + 1
                rtmp.append(rbltmp)
    r = sum(rtmp)
    print(r)
    return r

# (index, b, l, group, weight, importDay, exportDay)
def myRules(initDock, boxStream):
    num = len(boxStream) # 4494
    indexs = boxStream.index.values
    groups = boxStream.group.values
    weights = boxStream.weight.values
    importDays = boxStream.importDay.values
    exportDays = boxStream.exportDay.values
    # 按次序同组箱依次堆放
    n = 1 # 区段编号
    dock = initDock[n]
    br = list()
    for index in range(num):
        dbtmp = list()
        for nt in range(DN):
            dockt = initDock[nt + 1]
            for b in range(BN):
                isNotFull = False
                for l in range(LN):
                    if len(dockt[b][l]) < TN:
                        isNotFull = True
                        break
                    else:
                        continue
                groupn = groups[index]
                if isNotFull:
                    dgtmp = list()
                    for nn in range(DN):
                        dockn = initDock[nn + 1]
                        for bn in range(BN):
                            for ln in range(LN):
                                boxTien = dockn[bn][ln]
                                # 到所有同组箱的距离
                                for boxn in boxTien:
                                    if boxn[3] == groupn:
                                        dgtmp.append(abs(boxn[1] - b))
                    if dgtmp == []:
                        dgtmp.append(0)
                    dbtmp.append((nt, b, sum(dgtmp)))
        # 按距离排序升序
        dbtmp = sorted(dbtmp, key=lambda x: x[2])
        # 按贝位升序
        dbtmp = sorted(dbtmp, key=lambda x: x[1])
        # 按场区编号升序
        dbtmp = sorted(dbtmp, key=lambda x: x[0])
        docktb = initDock[dbtmp[0][0] + 1][dbtmp[0][1]]
        rbtmp = list()
        for l in range(LN):
            if len(docktb[l]) >= TN:
                continue
            else:
                rbltmp = 0
                boxTiett = docktb[l]
                bexportDays = list()
                bindexs = list()
                for iindex, box in enumerate(boxTiett):
                    bexportDays.append(box[6])
                    bindexs.append(iindex)
                    for bexportDay, bindex in zip(bexportDays, bindexs):
                        if box[6] > bexportDay:
                            rbltmp = rbltmp + iindex - bindex
                            break
                rbtmp.append((l, rbltmp))
                if rbtmp == []:
                    rbtmp.append((l, 0))
        rbtmp = sorted(rbtmp, key=lambda x: x[1])
        # print(dbtmp[0], rbtmp[0])
        br.append(rbtmp[0][1])
        box = (indexs[index], dbtmp[0][1], rbtmp[0][0], groups[index], weights[index], importDays[index], exportDays[index])
        dock = initDock[dbtmp[0][0] + 1]
        dock[dbtmp[0][1]][rbtmp[0][0]].append(box)
        dock = initDock[n]
        # 满
        if (index + 1) % (BN * LN * TN) == 0:
            n = -1 if n + 1 > DN else (n + 1)
            if n == -1:
                print('溢出')
                break
            dock = initDock[n]
    print('breshuffle:', sum(br))
    print(n)
    return distance(initDock=initDock), reshuffle(initDock=initDock)

# (index, i, j, group, weight, importDay, exportDay)
def rules(initDock, boxStream):
    num = len(boxStream) # 4494
    indexs = boxStream.index.values
    groups = boxStream.group.values
    weights = boxStream.weight.values
    importDays = boxStream.importDay.values
    exportDays = boxStream.exportDay.values
    # 按次序同组箱依次堆放
    n = 1 # 区段编号
    dock = initDock[n]
    for index in range(num):
        isAssig = False
        for b in range(BN):
            for l in range(LN):
                if len(dock[b][l]) >= TN: # * theta  
                    continue
                else:
                    isAssig = True
                    box = (indexs[index], b, l, groups[index], weights[index], importDays[index], exportDays[index])
                    dock[b][l].append(box)
                    break
            if isAssig:
                break
        # 满
        if (index + 1) % (BN * LN * TN) == 0:
            n = -1 if n + 1 > DN else (n + 1)
            if n == -1:
                print('溢出')
                break
            dock = initDock[n]
    # print(n)
    return distance(initDock=initDock), reshuffle(initDock=initDock)

if __name__ == "__main__":
    fakeData = pd.read_csv('fakeData.csv')
    initData(fakeData=fakeData)
    fakeData = pd.read_csv('fakeData.csv')
    print(fakeData)
    fakeB10 = fakeData[fakeData.exportDay <= 10]  # 前10天数据
    fakeA10 = fakeData[fakeData.exportDay > 10]  # 后20天数据
    # 按组内箱数量降序，按进场时间升序
    fakeB10 = pd.merge(fakeB10, fakeB10.groupby('group').size().reset_index(name='groupSize'), on='group')
    fakeB10 = fakeB10.sort_values(by=['groupSize'], ascending=False, kind='mergesort')
    fakeB10 = fakeB10.sort_values(by=['importDay'],  kind='mergesort')
    start = time.time()
    d, r = rules(initDock=initDock(fakeB10), boxStream=fakeB10)
    end = time.time()
    print(f'rules cost {end - start}.')
    start = time.time()
    md, mr = myRules(initDock=initDock(fakeB10), boxStream=fakeB10)
    end = time.time()
    print(f'myRules cost {end - start}.')
    f = open('log.txt', 'w+')
    f.write(f'r_distance:{d}, r_reshuffle:{r}\n')
    f.write(f'mr_distance:{md}, mr_reshuffle:{mr}\n')
    print('distance:', d, 'reshuffle:', r)
    print('distance:', md, 'reshuffle:', mr)
    print(f'd_improve:{(md - d) / d}, r_improve:{(mr - r) / r}')

    # # 初始化堆场
    # initD = initDock(boxStream=fakeB10)
    # # 按组内箱数量降序，按进场时间升序
    # fakeA10 = pd.merge(fakeA10, fakeA10.groupby('group').size().reset_index(name='groupSize'), on='group')
    # fakeA10 = fakeA10.sort_values(by=['groupSize'], ascending=False, kind='mergesort')
    # fakeA10 = fakeA10.sort_values(by=['importDay'],  kind='mergesort')
    # # 不初始化
    # # fakeData = pd.merge(fakeData, fakeData.groupby('group').size().reset_index(name='groupSize'), on='group')
    # # fakeData = fakeData.sort_values(by=['groupSize'], ascending=False, kind='mergesort')
    # # fakeData = fakeData.sort_values(by=['importDay'],  kind='mergesort')
    # start = time.time()
    # d, r = rules(initDock=initDock(fakeData), boxStream=fakeData)
    # end = time.time()
    # print(f'rules cost {end - start}.')
    # start = time.time()
    # md, mr = myRules(initDock=initDock(fakeData), boxStream=fakeData)
    # end = time.time()
    # print(f'myRules cost {end - start}.')
    # f = open('log.txt', 'w+')
    # f.write(f'r_distance:{d}, r_reshuffle:{r}\n')
    # f.write(f'mr_distance:{md}, mr_reshuffle:{mr}\n')
    # print('distance:', d, 'reshuffle:', r)
    # print('distance:', md, 'reshuffle:', mr)
    # print(f'd_improve:{(d - md) / d}, r_improve:{(r - mr) / r}')
