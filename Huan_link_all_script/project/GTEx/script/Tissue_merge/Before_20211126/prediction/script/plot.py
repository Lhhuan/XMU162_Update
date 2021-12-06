def Find_Optimal_Cutoff(TPR, FPR, threshold):
    y = TPR - FPR
    Youden_index = np.argmax(y)  # Only the first occurrence is returned.
    optimal_threshold = threshold[Youden_index]
    point = [FPR[Youden_index], TPR[Youden_index]]
    return optimal_threshold, point

def cm(prediction,ground_truth):
    prediction = prediction.copy()
    ground_truth = ground_truth.copy()
    ground_truth = ground_truth.astype(np.int32)
    if ground_truth.ndim == 2:
        ground_truth = ground_truth.squeeze()
    if prediction.shape[1] == 1:
        prediction_ = np.array(th.sigmoid(th.tensor(prediction)))
        fpr, tpr, thresholds = roc_curve(ground_truth, np.array(prediction_), )
        optimal_th, optimal_point = Find_Optimal_Cutoff(TPR=tpr, FPR=fpr, threshold=thresholds)
        prediction[prediction_>=optimal_th]=1
        prediction[prediction_ < optimal_th] = 0
        prediction = prediction.squeeze().astype(np.int32)
    else:
        print('done')
        prediction_ = np.array(th.softmax(th.tensor(prediction),axis=1)[:,1])
        fpr, tpr, thresholds = roc_curve(ground_truth, prediction_, )
        optimal_th, optimal_point = Find_Optimal_Cutoff(TPR=tpr, FPR=fpr, threshold=thresholds)
        prediction_[prediction_>=optimal_th]=1
        prediction_[prediction_ < optimal_th] = 0
        prediction = prediction_.squeeze().astype(np.int32)
        print(prediction.shape)
        # prediction = np.argmax(np.array(prediction),axis=1)
    x =~(ground_truth ^ prediction) + 2 #1为预测对，0为错

    print(np.sum(ground_truth==1),np.sum(ground_truth==0))
    TP = np.sum((ground_truth==1) & (x==1))/ np.sum(ground_truth==1)
    TN = np.sum((ground_truth==0) & (x==1))/ np.sum(ground_truth==0)
    FP = np.sum((ground_truth==0) & (x==0))/ np.sum(ground_truth==0)
    FN = np.sum((ground_truth==1) & (x==0))/ np.sum(ground_truth==1)

    accuracy = (TP+TN) / (TP+TN+FP+FN)
    precision = TP/(TP+FP)
    recall = TP/(TP+FN) #sensitivity
    specificity = TN/(TN+FP)


    F1_score = 2*(precision*recall)/(precision+recall)