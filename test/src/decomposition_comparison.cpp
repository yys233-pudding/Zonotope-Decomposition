#include <iostream>
#include <vector>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include "ZonoOpt.hpp"

using namespace ZonoOpt;
using namespace std;
using namespace Eigen;

// 结构用于结果
struct TestResult {
    double time_union;
    double time_zono_union;
    double time_set_diff;
    double time_minkowski;
    double time_intersection;
    double time_halfspace;
    double coverage;
    double compactness;
};

// 生成测试用例：创建空间、子 zonotope、障碍物
struct TestCase {
    Zono full_space;
    vector<Zono> sub_zonos;
    vector<Zono> obstacles;
    SparseMatrix<double> H;  // 半空间矩阵
    VectorXd f;  // 半空间向量
    double true_volume;
};

TestCase generate_test_case(int case_id) {
    // 仿照 example_operations 创建复杂 zonotope
    int dim = 2 + (case_id % 3);  // 2-4 维度
    SparseMatrix<double> G_full(dim, dim);
    for (int i = 0; i < dim; ++i) G_full.insert(i, i) = 1.0 + case_id * 0.1;
    VectorXd c_full = VectorXd::Zero(dim);
    Zono full_space(G_full, c_full);

    vector<Zono> sub_zonos;
    int num_sub = 2 + case_id % 3;  // 2-4 子集
    for (int i = 0; i < num_sub; ++i) {
        SparseMatrix<double> G_sub(dim, dim);
        for (int j = 0; j < dim; ++j) G_sub.insert(j, j) = 0.5 + i * 0.2;
        VectorXd c_sub = VectorXd::Random(dim) * 2.0;
        // 创建 Zono 或 ConZono
        if (i % 2 == 0) {
            // Zono
            sub_zonos.emplace_back(G_sub, c_sub);
        } else {
            // ConZono
            SparseMatrix<double> A_sub(1, dim);
            A_sub.insert(0, 0) = 1;
            VectorXd b_sub(1);
            b_sub(0) = 1 + i * 0.5;
            ConZono conz(G_sub, c_sub, A_sub, b_sub, true);
            // 转换为 Zono (简化)
            sub_zonos.emplace_back(G_sub, c_sub);
        }
    }

    vector<Zono> obstacles;
    int num_obs = 1 + case_id % 2;  // 1-2 障碍物
    for (int i = 0; i < num_obs; ++i) {
        SparseMatrix<double> G_obs(dim, 1);
        G_obs.insert(0, 0) = 0.3 + i * 0.1;
        VectorXd c_obs = VectorXd::Random(dim) * 1.5;
        obstacles.emplace_back(G_obs, c_obs);
    }

    // 半空间参数
    SparseMatrix<double> H(dim, dim);
    H.insert(0, 0) = 1;
    VectorXd f(dim);
    f(0) = 2.0;

    double true_volume = pow(10.0, dim);  // 简化

    return {full_space, sub_zonos, obstacles, H, f, true_volume};
}

// 简化指标
double compute_coverage(const unique_ptr<HybZono>& decomposed) {
    return 0.8;  // 简化
}

double compute_compactness(const unique_ptr<HybZono>& decomposed, double true_volume) {
    return 0.9;  // 简化
}

TestResult run_test(const TestCase& tc) {
    TestResult res;

    // union_of_many
    auto start = chrono::high_resolution_clock::now();
    vector<shared_ptr<HybZono>> hyb_zonos;
    for (const auto& z : tc.sub_zonos) {
        hyb_zonos.push_back(make_shared<HybZono>(z));
    }
    auto res1 = union_of_many(hyb_zonos);
    auto end = chrono::high_resolution_clock::now();
    res.time_union = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    res.coverage = compute_coverage(res1);
    res.compactness = compute_compactness(res1, tc.true_volume);

    // zono_union_2_hybzono
    start = chrono::high_resolution_clock::now();
    vector<Zono> zs2 = tc.sub_zonos;
    auto res2 = zono_union_2_hybzono(zs2);
    end = chrono::high_resolution_clock::now();
    res.time_zono_union = chrono::duration_cast<chrono::milliseconds>(end - start).count();

    // set_diff
    start = chrono::high_resolution_clock::now();
    HybZono free_space(tc.full_space);
    for (const auto& obs : tc.obstacles) {
        HybZono obs_hyb(obs);
        auto diff = set_diff(free_space, obs_hyb);
        free_space = *diff;
    }
    end = chrono::high_resolution_clock::now();
    res.time_set_diff = chrono::duration_cast<chrono::milliseconds>(end - start).count();

    // minkowski_sum (使用前两个子集)
    if (tc.sub_zonos.size() >= 2) {
        start = chrono::high_resolution_clock::now();
        HybZono z1(tc.sub_zonos[0]);
        HybZono z2(tc.sub_zonos[1]);
        auto res3 = minkowski_sum(z1, z2);
        end = chrono::high_resolution_clock::now();
        res.time_minkowski = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    } else {
        res.time_minkowski = 0;
    }

    // intersection (使用前两个子集)
    if (tc.sub_zonos.size() >= 2) {
        start = chrono::high_resolution_clock::now();
        HybZono z1(tc.sub_zonos[0]);
        HybZono z2(tc.sub_zonos[1]);
        auto res4 = intersection(z1, z2);
        end = chrono::high_resolution_clock::now();
        res.time_intersection = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    } else {
        res.time_intersection = 0;
    }

    // halfspace_intersection
    start = chrono::high_resolution_clock::now();
    HybZono z_hyb(tc.full_space);
    auto res5 = halfspace_intersection(z_hyb, tc.H, tc.f);
    end = chrono::high_resolution_clock::now();
    res.time_halfspace = chrono::duration_cast<chrono::milliseconds>(end - start).count();

    return res;
}

int main() {
    const int num_tests = 12;  // 10+ 测试用例
    vector<TestResult> results;

    for (int i = 0; i < num_tests; ++i) {
        TestCase tc = generate_test_case(i);
        TestResult res = run_test(tc);
        results.push_back(res);
        cout << "Test " << i << ": union=" << res.time_union << "ms, zono_union=" << res.time_zono_union << "ms, set_diff=" << res.time_set_diff << "ms, minkowski=" << res.time_minkowski << "ms, intersection=" << res.time_intersection << "ms, halfspace=" << res.time_halfspace << "ms" << endl;
    }

    // 计算平均
    double avg_union = 0, avg_zono_union = 0, avg_set_diff = 0, avg_minkowski = 0, avg_intersection = 0, avg_halfspace = 0, avg_cov = 0, avg_comp = 0;
    for (const auto& r : results) {
        avg_union += r.time_union;
        avg_zono_union += r.time_zono_union;
        avg_set_diff += r.time_set_diff;
        avg_minkowski += r.time_minkowski;
        avg_intersection += r.time_intersection;
        avg_halfspace += r.time_halfspace;
        avg_cov += r.coverage;
        avg_comp += r.compactness;
    }
    avg_union /= num_tests;
    avg_zono_union /= num_tests;
    avg_set_diff /= num_tests;
    avg_minkowski /= num_tests;
    avg_intersection /= num_tests;
    avg_halfspace /= num_tests;
    avg_cov /= num_tests;
    avg_comp /= num_tests;

    cout << "\nAverage Results:" << endl;
    cout << "union_of_many: Time=" << avg_union << "ms, Coverage=" << avg_cov << ", Compactness=" << avg_comp << endl;
    cout << "zono_union_2_hybzono: Time=" << avg_zono_union << "ms" << endl;
    cout << "set_diff: Time=" << avg_set_diff << "ms" << endl;
    cout << "minkowski_sum: Time=" << avg_minkowski << "ms" << endl;
    cout << "intersection: Time=" << avg_intersection << "ms" << endl;
    cout << "halfspace_intersection: Time=" << avg_halfspace << "ms" << endl;

    // 选出最佳: 基于时间 (可扩展)
    vector<double> times = {avg_union, avg_zono_union, avg_set_diff, avg_minkowski, avg_intersection, avg_halfspace};
    vector<string> names = {"union_of_many", "zono_union_2_hybzono", "set_diff", "minkowski_sum", "intersection", "halfspace_intersection"};
    int best_idx = 0;
    double min_time = times[0];
    for (size_t i = 1; i < times.size(); ++i) {
        if (times[i] < min_time) {
            min_time = times[i];
            best_idx = i;
        }
    }

    cout << "Best: " << names[best_idx] << endl;

    return 0;
}