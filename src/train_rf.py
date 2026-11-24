import numpy as np
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error

X_train, y_train = np.load('data/X_train.npy'), np.load('data/y_train.npy')
X_test, y_test = np.load('data/X_test.npy'), np.load('data/y_test.npy')

def evaluate(name, model):
    model.fit(X_train, y_train)          # öğrenme: sadece train verisini görür
    pred = model.predict(X_test)         # tahmin: hiç görmediği test molekülleri
    r2 = r2_score(y_test, pred)
    rmse = np.sqrt(mean_squared_error(y_test, pred))
    print(f'{name:15s} R2: {r2:.3f}  RMSE: {rmse:.3f}')
    return pred

evaluate('Baseline (ort.)', DummyRegressor(strategy='mean'))

rf = RandomForestRegressor(n_estimators=300, n_jobs=-1, random_state=42)
pred_rf = evaluate('Random Forest', rf)
np.save('data/pred_rf.npy', pred_rf)