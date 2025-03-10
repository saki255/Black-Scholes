#include <iostream>
#include <cmath>
#include <iomanip>

// 数値計算のための定数
const double PI = 3.14159265358979323846;

// 標準正規分布の累積分布関数
double normalCDF(double x) {
    // 近似計算
    if (std::abs(x) > 6.0) {
        return x < 0.0 ? 0.0 : 1.0; // 極端な値の処理
    }
    
    double k = 1.0 / (1.0 + 0.2316419 * std::abs(x));
    double z = 1.0 / std::sqrt(2.0 * PI) * std::exp(-0.5 * x * x);
    double a = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));
    
    if (x >= 0.0) {
        return 1.0 - z * a;
    } else {
        return z * a;
    }
}

// 日経225用ブラック・ショールズモデルによるコールオプション価格計算
double nikkeiCallOption(double S, double K, double r, double dividendYield, double sigma, double T) {
    // S: 日経平均株価
    // K: 権利行使価格
    // r: 無リスク金利
    // dividendYield: 配当利回り
    // sigma: ボラティリティ
    // T: 満期までの期間（年）
    
    double d1 = (std::log(S / K) + (r - dividendYield + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    
    return S * std::exp(-dividendYield * T) * normalCDF(d1) - K * std::exp(-r * T) * normalCDF(d2);
}

// 日経225用ブラック・ショールズモデルによるプットオプション価格計算
double nikkeiPutOption(double S, double K, double r, double dividendYield, double sigma, double T) {
    double d1 = (std::log(S / K) + (r - dividendYield + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    
    return K * std::exp(-r * T) * normalCDF(-d2) - S * std::exp(-dividendYield * T) * normalCDF(-d1);
}

// コールオプションのデルタ計算
double callDelta(double S, double K, double r, double dividendYield, double sigma, double T) {
    double d1 = (std::log(S / K) + (r - dividendYield + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return std::exp(-dividendYield * T) * normalCDF(d1);
}

// プットオプションのデルタ計算
double putDelta(double S, double K, double r, double dividendYield, double sigma, double T) {
    double d1 = (std::log(S / K) + (r - dividendYield + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    return std::exp(-dividendYield * T) * (normalCDF(d1) - 1.0);
}

// ガンマ計算（コール・プット共通）
double gamma(double S, double K, double r, double dividendYield, double sigma, double T) {
    double d1 = (std::log(S / K) + (r - dividendYield + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double pdf = (1.0 / std::sqrt(2.0 * PI)) * std::exp(-0.5 * d1 * d1); // 標準正規分布の確率密度関数
    
    return std::exp(-dividendYield * T) * pdf / (S * sigma * std::sqrt(T));
}

// ベガ計算（コール・プット共通、1%変化あたり）
double vega(double S, double K, double r, double dividendYield, double sigma, double T) {
    double d1 = (std::log(S / K) + (r - dividendYield + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double pdf = (1.0 / std::sqrt(2.0 * PI)) * std::exp(-0.5 * d1 * d1);
    
    return S * std::exp(-dividendYield * T) * std::sqrt(T) * pdf / 100.0; // 1%変化あたりに調整
}

// シータ計算（日次）
double callTheta(double S, double K, double r, double dividendYield, double sigma, double T) {
    double d1 = (std::log(S / K) + (r - dividendYield + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    double pdf = (1.0 / std::sqrt(2.0 * PI)) * std::exp(-0.5 * d1 * d1);
    
    double theta = -S * std::exp(-dividendYield * T) * pdf * sigma / (2 * std::sqrt(T))
                  + dividendYield * S * std::exp(-dividendYield * T) * normalCDF(d1)
                  - r * K * std::exp(-r * T) * normalCDF(d2);
    
    return theta / 365.0; // 年間から日次に変換
}

double putTheta(double S, double K, double r, double dividendYield, double sigma, double T) {
    double d1 = (std::log(S / K) + (r - dividendYield + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);
    double pdf = (1.0 / std::sqrt(2.0 * PI)) * std::exp(-0.5 * d1 * d1);
    
    double theta = -S * std::exp(-dividendYield * T) * pdf * sigma / (2 * std::sqrt(T))
                  - dividendYield * S * std::exp(-dividendYield * T) * normalCDF(-d1)
                  + r * K * std::exp(-r * T) * normalCDF(-d2);
    
    return theta / 365.0; // 年間から日次に変換
}

int main() {
    // 日経225オプション用のパラメータ設定
    double S = 36420.0;         // 日経平均株価（例: 25,000円）
    double K = 38000.0;         // 権利行使価格
    double r = 0.0478;           // 無リスク金利（例: 0.5%）
    double dividendYield = 0.02; // 配当利回り（例: 1.5%）
    double sigma = 0.23;        // ボラティリティ（例: 20%）
    double T = 7.0 / 365.0;    // 満期までの期間（例: 30日）
    
    // オプション価格計算
    double call_price = nikkeiCallOption(S, K, r, dividendYield, sigma, T);
    double put_price = nikkeiPutOption(S, K, r, dividendYield, sigma, T);
    
    // グリークス計算
    double call_delta = callDelta(S, K, r, dividendYield, sigma, T);
    double put_delta = putDelta(S, K, r, dividendYield, sigma, T);
    double option_gamma = gamma(S, K, r, dividendYield, sigma, T);
    double option_vega = vega(S, K, r, dividendYield, sigma, T);
    double call_theta = callTheta(S, K, r, dividendYield, sigma, T);
    double put_theta = putTheta(S, K, r, dividendYield, sigma, T);
    
    // 結果の表示
    std::cout << "日経225オプション価格計算（ブラック・ショールズモデル）" << std::endl;
    std::cout << "=================================================" << std::endl;
    std::cout << "パラメータ:" << std::endl;
    std::cout << "日経平均株価: " << S << " 円" << std::endl;
    std::cout << "権利行使価格: " << K << " 円" << std::endl;
    std::cout << "無リスク金利: " << r * 100 << " %" << std::endl;
    std::cout << "配当利回り: " << dividendYield * 100 << " %" << std::endl;
    std::cout << "ボラティリティ: " << sigma * 100 << " %" << std::endl;
    std::cout << "満期までの期間: " << T * 365 << " 日" << std::endl;
    std::cout << "=================================================" << std::endl;
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "オプション価格:" << std::endl;
    std::cout << "コールオプション: " << call_price << " 円" << std::endl;
    std::cout << "プットオプション: " << put_price << " 円" << std::endl;
    std::cout << "=================================================" << std::endl;
    
    std::cout << "グリークス:" << std::endl;
    std::cout << "コールデルタ: " << call_delta << std::endl;
    std::cout << "プットデルタ: " << put_delta << std::endl;
    std::cout << "ガンマ: " << option_gamma << std::endl;
    std::cout << "ベガ (1%変化あたり): " << option_vega << " 円" << std::endl;
    std::cout << "コールシータ (1日): " << call_theta << " 円" << std::endl;
    std::cout << "プットシータ (1日): " << put_theta << " 円" << std::endl;
    std::cout << "=================================================" << std::endl;
    
    // プット・コール・パリティの確認
    double parity_check = call_price - put_price - S * std::exp(-dividendYield * T) + K * std::exp(-r * T);
    std::cout << "プット・コール・パリティ誤差: " << std::abs(parity_check) << " 円" << std::endl;
    
    // ユーザー入力を待機
    std::cout << "\nEnterキーを押して終了...";
    std::cin.get();
    
    return 0;
}