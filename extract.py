
import sys

import numpy as np

filename = sys.argv[1]   # tig00000007
lines = int(sys.argv[2])

ratioThresh = int(sys.argv[3])    # 2
depthThresh = int(sys.argv[4])   # 30
posThresh = int(sys.argv[5])     # 63




# filename = 'tig00000331.peak'
# lines = 145154
#
# ratioThresh = 2
# depthThresh = 15
# posThresh = 63


def findPeak(filename, lines, ratioThresh=2, posThresh=63):
    with open(filename) as f:
        result = []
        record = []
        tempRecord = []
        name = filename
        rise = []
        drop = []


        for i, line in enumerate(f):

            # process first 70 lines
            if i < 70:
                name, currPos, _, currDepth = line.split()
                tempRecord.append((int(currPos), int(currDepth)))
                # add to record
                if i == 69:
                    prevDepth = np.mean(np.asarray(tempRecord)[0:7, 1])

                    currDepth = np.mean(np.asarray(tempRecord)[60:67, 1])

                    currPos = tempRecord[60][0]
                    print('currDepth= %s' % currDepth)
                    print('prevDepth= %s' % prevDepth)
                    print('line number = %d' % i)

                    rise.append((currPos,prevDepth,currDepth))
                    record.append((0,currPos,int(prevDepth),int(currDepth)))

                    # add startPos and depth to record at the beginning
                    #record.append((currPos, currDepth))
                    print('rise append pos=%d prevDepth= %d currDepth= %d ' % (currPos, prevDepth,currDepth))
                    print('First')
                    #first = True
                continue

            # process last line
            if i == lines-70:
                # if record:
                #     print('currDepth= %s' % currDepth)
                #     print('prevDepth= %s' % prevDepth)
                #     print('line number = %d' % i)
                #     print('startPos = %d , endPos = %d' % (startPos,currPos))
                #     _, endPos, _ = line.split()
                #     startPos, _ = record[-1]
                #     result.append((startPos, endPos))
                resultFile = str(filename + '.peak')
                f2 = open(resultFile, 'w')
                for mark,pos,pre,cur in record:
                    f2.write(str(mark) + ' ' + str(pos) + ' ' + str(pre) + ' ' + str(cur) + '\n')
                return

            # Process normal lines
            _, pos, _, depth = line.split()
            tempRecord.pop(0)
            tempRecord.append((int(pos), int(depth)))

            prevDepth = np.mean(np.asarray(tempRecord)[0:7, 1])

            currDepth = np.mean(np.asarray(tempRecord)[60:67, 1])

            prevStd = np.std(np.asarray(tempRecord)[0:7, 1])

            currStd = np.std(np.asarray(tempRecord)[60:67, 1])

            currPos = tempRecord[60][0]

            ratio = currDepth / prevDepth

            # when coverage rise
            if ratio > ratioThresh and prevStd < 15:

                if rise:
                    prevPos = rise[-1][0]
                    # rec_prevDepth = rise[-1][1]
                    if abs(prevPos - currPos) > posThresh:
                        print('coverage rise')
                        print('ratio= %f' % ratio)
                        print('currDepth= %s' % currDepth)
                        print('prevDepth= %s' % prevDepth)
                        print('line number = %d' % i)
                        print('rise append pos=%d prevDepth= %d currDepth= %d' % (currPos, prevDepth,currDepth))
                        rise.append((currPos, prevDepth,currDepth))
                        record.append((0, currPos, int(prevDepth), int(currDepth)))
                else:
                    print('coverage rise')
                    print('ratio= %f' % ratio)
                    print('currDepth= %s' % currDepth)
                    print('prevDepth= %s' % prevDepth)
                    print('line number = %d' % i)
                    print('rise append pos=%d prevDepth= %d currDepth= %d' % (currPos, prevDepth, currDepth))
                    rise.append((currPos, prevDepth, currDepth))
                    record.append((0, currPos, int(prevDepth), int(currDepth)))

                continue

            # when coverage drop
            if ratio < 1 / ratioThresh and currStd < 15:
                if drop:
                    prevPos = drop[-1][0]
                    if abs(prevPos - currPos) > posThresh:
                        print('coverage drop')
                        print('ratio= %f' % ratio)
                        print('currDepth= %s' % currDepth)
                        print('prevDepth= %s' % prevDepth)
                        print('line number = %d' % i)
                        print('drop append pos=%d prevDepth= %d currDepth= %d' % (currPos, prevDepth, currDepth))
                        drop.append((currPos, prevDepth, currDepth))
                        record.append((1, currPos, int(prevDepth), int(currDepth)))

                else:
                    print('coverage drop')
                    print('ratio= %f' % ratio)
                    print('currDepth= %s' % currDepth)
                    print('prevDepth= %s' % prevDepth)
                    print('line number = %d' % i)
                    print('drop append pos=%d prevDepth= %d currDepth= %d' % (currPos, prevDepth, currDepth))
                    drop.append((currPos, prevDepth, currDepth))
                    record.append((1, currPos, int(prevDepth), int(currDepth)))


def mergePeak(filename, depthThresh = 30):

    with open(str(filename+'.peak')) as f:

        data = np.fromfile(f, dtype=int, sep=' ')
        data = np.reshape(data, (-1, 4))

        lastSig = 0
        lastPos = 0
        result = []
        lines = {}
        skip = False

        for i in range(data.shape[0]):
            sig, pos, preDepth, curDepth= data[i]

            if i == data.shape[0]-1:
                if len(result) == 0:
                    return
                with open(str(filename + '.bed'),'w') as f1:
                    for j in range(1,len(result),2):
                        f1.write(filename + ' ' + str(result[j-1][1]) + ' ' + str(result[j][1]) + '\n')
                    return
            if preDepth < depthThresh or curDepth < depthThresh:
                continue
            lines[pos] = i
            if i == 0 and sig == 0:
                lastPos = pos
                lastSig = sig

            change = (sig != lastSig)
            if change:
                if sig == 1:
                    if abs(pos - lastPos) < 10000:
                        if lastPos in lines:
                            result.append(data[lines[lastPos]])
                    else:
                        skip = True      # too long repeat, maybe no match rise

                else:
                    if not skip:  # have rise before
                        if lastPos in lines:
                            result.append(data[lines[lastPos]]) # append last drop position

                    skip = False

                lastSig = sig
                lastPos = pos

            else:

                if sig == 0:
                    if abs(pos - lastPos) > 10000:    # two consecutive rise with long space, update to the second one

                        lastPos = pos
                        lastSig = sig
                    else:
                        if lastPos == 0:
                            lastPos = pos
                            lastSig = sig

                else:
                    if not skip:    # have rise before
                        if abs(pos - lastPos) < 10000:
                            lastPos = pos
                            lastSig = sig


# findPeak(filename,lines,ratioThresh,posThresh)
# mergePeak(filename,depthThresh)






def findPeak_new(filename, lines, ratioThresh=2, posThresh=63):
    with open(filename) as f:

        record = []

        rise = []
        drop = []

        data = np.fromfile(f, dtype=int, sep=' ')

        data = np.reshape(data, (-1, 3))


        # init

        prevDepth = np.mean(data[0:7, 2])

        currDepth = np.mean(data[63:70, 2])

        currPos = data[63,1]
        print('currDepth= %s' % currDepth)
        print('prevDepth= %s' % prevDepth)
        print('First')


        rise.append((currPos, prevDepth, currDepth))
        record.append((0, currPos, int(prevDepth), int(currDepth)))

        for i in range(71,data.shape[0]-70,2):

            prevDepth = np.mean(data[i-70:i-63, 2])

            currDepth = np.mean(data[i-7:i, 2])

            prevStd = np.std(data[i-70:i-63, 2])

            currStd = np.std(data[i-7:i, 2])

            currPos = data[i-7,0]

            ratio = currDepth / prevDepth

            # when coverage rise
            if ratio > ratioThresh and prevStd < 15:

                if rise:
                    prevPos = rise[-1][0]
                    # rec_prevDepth = rise[-1][1]
                    if abs(prevPos - currPos) > posThresh:
                        print('coverage rise')
                        print('ratio= %f' % ratio)
                        print('currDepth= %s' % currDepth)
                        print('prevDepth= %s' % prevDepth)
                        print('line number = %d' % i)
                        print('rise append pos=%d prevDepth= %d currDepth= %d' % (currPos, prevDepth,currDepth))
                        rise.append((currPos, prevDepth,currDepth))
                        record.append((0, currPos, int(prevDepth), int(currDepth)))
                else:
                    print('coverage rise')
                    print('ratio= %f' % ratio)
                    print('currDepth= %s' % currDepth)
                    print('prevDepth= %s' % prevDepth)
                    print('line number = %d' % i)
                    print('rise append pos=%d prevDepth= %d currDepth= %d' % (currPos, prevDepth, currDepth))
                    rise.append((currPos, prevDepth, currDepth))
                    record.append((0, currPos, int(prevDepth), int(currDepth)))

                continue

            # when coverage drop
            if ratio < 1 / ratioThresh and currStd < 15:
                if drop:
                    prevPos = drop[-1][0]
                    if abs(prevPos - currPos) > posThresh:
                        print('coverage drop')
                        print('ratio= %f' % ratio)
                        print('currDepth= %s' % currDepth)
                        print('prevDepth= %s' % prevDepth)
                        print('line number = %d' % i)
                        print('drop append pos=%d prevDepth= %d currDepth= %d' % (currPos, prevDepth, currDepth))
                        drop.append((currPos, prevDepth, currDepth))
                        record.append((1, currPos, int(prevDepth), int(currDepth)))

                else:
                    print('coverage drop')
                    print('ratio= %f' % ratio)
                    print('currDepth= %s' % currDepth)
                    print('prevDepth= %s' % prevDepth)
                    print('line number = %d' % i)
                    print('drop append pos=%d prevDepth= %d currDepth= %d' % (currPos, prevDepth, currDepth))
                    drop.append((currPos, prevDepth, currDepth))
                    record.append((1, currPos, int(prevDepth), int(currDepth)))

        resultFile = str(filename + '.peak')
        f2 = open(resultFile, 'w')

        for mark, pos, pre, cur in record:
            f2.write(str(mark) + ' ' + str(pos) + ' ' + str(pre) + ' ' + str(cur) + '\n')

        return

findPeak_new(filename,lines,ratioThresh,posThresh)
mergePeak(filename,depthThresh)