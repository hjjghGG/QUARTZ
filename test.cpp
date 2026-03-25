#include <climits>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <iostream>
#include <string>
#include <regex>
#include <sys/timeb.h>


#include "KC_pu.h"


using namespace std;


#ifndef CLK_PER_SEC
#ifdef CLOCKS_PER_SEC
#define CLK_PER_SEC CLOCKS_PER_SEC
#endif
#endif


double _theta = 0.01;
double alpha = 0.01;
int max_buckets = 2048;
#define LOOP_NUM 1
#define QUERY_PER_UPDATE 10000000
#define QUERY_PER_UPDATESMALL 10000


int CalculateMaxBucketsForMemory(int total_mem, int k) {
    if (total_mem <= 0 || k <= 0) {
        std::cerr << "Error: total_mem and k must be positive values" << std::endl;
        return -1;
    }


    const int SIZEOF_KCU_TYPE = sizeof(KCU_type);
    const int SIZEOF_KCU_ITEM = sizeof(KCU_item);
    const int SIZEOF_KCU_ITEM_PTR = sizeof(KCU_item *);
    const int UDD_PER_K = (56 + 56 + 32); // psketch + mapstore + mapping

    int kcu_fixed = SIZEOF_KCU_TYPE;

    int kcu_per_k = SIZEOF_KCU_ITEM + SIZEOF_KCU_ITEM_PTR +
                    (4 + 1) * KCU_CELL_SIZE; // (unsigned int + uint8_t) * KCU_CELL_SIZE

    int udd_part = UDD_PER_K * k;

    int fixed_cost = kcu_fixed + (kcu_per_k * k) + udd_part;

    if (fixed_cost >= total_mem) {
        std::cerr << "Error: Fixed cost (" << fixed_cost
                << " bytes) already exceeds total memory ("
                << total_mem << " bytes)" << std::endl;
        return 0;
    }

    int available_for_buckets = total_mem - fixed_cost;
    const int BUCKET_SIZE = 12;

    int max_buckets = available_for_buckets / BUCKET_SIZE;

    std::cout << "  Calculated max_buckets: " << max_buckets << std::endl;

    return max_buckets;
}

double computePercentile(const std::vector<double> &vec, double item) {
    double quantile = 0.0;

    if (item <= vec.front()) {
        quantile = 0.0;
    } else if (item >= vec.back()) {
        quantile = 1.0;
    } else {
        auto it = std::lower_bound(vec.begin(), vec.end(), item);
        size_t pos = it - vec.begin();

        if (pos < vec.size() && vec[pos] == item) {
            while (pos > 0 && vec[pos - 1] == item) {
                --pos;
            }
            quantile = static_cast<double>(pos) / vec.size();
        } else {
            double dist_left = (pos == 0) ? INFINITY : item - vec[pos - 1];
            double dist_right = (pos == vec.size()) ? INFINITY : vec[pos] - item;

            if (dist_left <= dist_right) {
                quantile = static_cast<double>(pos - 1) / vec.size();
            } else {
                quantile = static_cast<double>(pos) / vec.size();
            }
        }
    }
    return quantile;
}

double query_for_val_vector(const vector<double>& vec, double rank) {
    vector<double> sorted_vec = vec;
    std::sort(sorted_vec.begin(), sorted_vec.end());
    int index = sorted_vec.size() * rank;

    if (index >= sorted_vec.size()) {
        index = sorted_vec.size() - 1;
    }

    return sorted_vec[index];
}

time_t now = time(0);

// Structure to store trace data
struct TraceData {
    int totalLines = 0;
    int traceSize = 0;
    int nonZerosLatency = 0;
    unordered_map<uint32_t, vector<double> > flowLatencies;
    unordered_map<uint32_t, uint32_t> flowCounters;
    unordered_map<uint32_t, uint32_t> lastUpdatedEpoch;
    vector<pair<uint32_t, double> > allPackets;
};

// Read trace data file
TraceData readTraceFile(const std::string &filename) {
    TraceData data;

    ifstream infile(filename);
    if (!infile.is_open()) {
        cout << "can't open file " << filename << endl;
        return data;
    }

    string line;
    while (getline(infile, line)) {
        if (line.empty()) continue;
        data.totalLines++;

        if (data.totalLines % 1000000 == 0) {
            cout << "processed " << data.totalLines << " lines" << endl;
        }

        std::vector<std::string> splits;
        std::string split;
        std::istringstream iss(line);
        char delim = ',';

        if (filename.find("case") != std::string::npos) delim = ' ';

        while (std::getline(iss, split, delim)) {
            splits.push_back(split);
        }

        if (splits.size() < 2) continue;

        try {
            uint32_t ID = (uint32_t) stoull(splits[0]);
            double latency = stoi(splits[1]);

            if (latency > 0) {
                data.nonZerosLatency++;
                data.allPackets.emplace_back(ID, latency);
                data.traceSize++;

                data.flowLatencies[ID].push_back(latency);

                int epoch = data.totalLines / KCU_EPOCH_SIZE + 1;
                auto it = data.lastUpdatedEpoch.find(ID);
                if (it == data.lastUpdatedEpoch.end() || it->second != epoch) {
                    data.flowCounters[ID]++;
                    data.lastUpdatedEpoch[ID] = epoch;
                }
            }
        } catch (const std::exception &e) {
            cout << "line: " << data.totalLines << " transform error " << e.what() << ", content: " << line << endl;
            continue;
        }
    }

    infile.close();
    cout << "Read " << data.totalLines << " lines, valid: " << data.traceSize
            << ", Non-zero: " << data.nonZerosLatency << endl;

    return data;
}

// Populate KCU sketch with pre-loaded data
void populateKCUSketch(KCU_type *KCU, const TraceData &data) {
    cout << "Populating KCU sketch with " << data.allPackets.size() << " packets..." << endl;

    for (const auto &packet: data.allPackets) {
        KCU_UpdateLatency(KCU, packet.first, packet.second);
    }
}

// Main test function using pre-loaded data
int Tail_latency_quantile_test_with_data(const std::string &filename, KCU_type *KCU_type_ptr, const TraceData &data) {
    cout << "debug: " << filename << endl;

    cout << "[DEBUG] Algorithm parameters:" << endl;
    cout << "  _theta = " << _theta << endl;
    cout << "  alpha = " << alpha << endl;
    cout << "  max_buckets = " << max_buckets << endl;
    cout << "  KCU_EPOCH_SIZE = " << KCU_EPOCH_SIZE << endl;
    cout << "  KCU_CELL_SIZE = " << KCU_CELL_SIZE << endl;

    ofstream outfile;
    outfile.open("../output/output.txt", ios::app);
    char *dt = ctime(&now);
    outfile << "==========================================" << endl;
    outfile << "Test info: " << dt;
    outfile << "parameter: theta=" << _theta << ", alpha=" << alpha << ", max_buckets=" << max_buckets << endl;
    outfile << "filename: " << filename << endl;
    outfile << "==========================================" << endl << endl;

    cout << "Result: " << data.totalLines << ", valid: " << data.traceSize
            << ", Non-zero: " << data.nonZerosLatency << endl;
    outfile << "Result: " << data.totalLines << ", valid: " << data.traceSize
            << ", Non-zero: " << data.nonZerosLatency << endl;

    cout << "\n[DEBUG] Flow statistics:" << endl;
    cout << "  Total unique flows: " << data.flowLatencies.size() << endl;

    vector<pair<uint32_t, size_t> > flow_freq;
    for (const auto &it: data.flowLatencies) {
        flow_freq.emplace_back(it.first, it.second.size());
    }
    sort(flow_freq.begin(), flow_freq.end(),
         [](const pair<uint32_t, size_t> &a, const pair<uint32_t, size_t> &b) {
             return a.second > b.second;
         });

    cout << "  Top 10 most frequent flows:" << endl;
    for (int i = 0; i < min(10, (int) flow_freq.size()); i++) {
        cout << "    Flow " << flow_freq[i].first << ": " << flow_freq[i].second << " packets" << endl;
    }

    cout << "\n[DEBUG] CountersMap statistics:" << endl;
    vector<pair<uint32_t, uint32_t> > epoch_counts;
    for (const auto &it: data.flowCounters) {
        epoch_counts.emplace_back(it.first, it.second);
    }
    sort(epoch_counts.begin(), epoch_counts.end(),
         [](const pair<uint32_t, uint32_t> &a, const pair<uint32_t, uint32_t> &b) {
             return a.second > b.second;
         });

    int total_epochs = data.totalLines / KCU_EPOCH_SIZE + 1;
    cout << "  Total epochs: " << total_epochs << endl;
    cout << "  Epoch threshold (theta * epochs): " << (_theta * total_epochs) << endl;
    cout << "  Top 10 flows by epoch count:" << endl;
    for (int i = 0; i < min(10, (int) epoch_counts.size()); i++) {
        bool above_threshold = epoch_counts[i].second > (_theta * total_epochs);
    }

    int thresh = ceil(_theta * data.traceSize);
    cout << " theta=" << _theta << ", traceSize=" << data.traceSize << ", thresh=" << thresh << endl;
    outfile << " theta=" << _theta << ", traceSize=" << data.traceSize << ", thresh=" << thresh << endl;

    std::pair<int, std::map<uint32_t, uint32_t> > result = KCU_Output2(KCU_type_ptr);
    std::map<uint32_t, uint32_t> bucketElements = result.second;
    cout << "ending_buckets: " << result.first << endl;
    cout << "\n[DEBUG] KCU_Output2 results:" << endl;
    cout << "  bucketElements size: " << bucketElements.size() << endl;

    cout << "Collision times: "
#ifdef ORG_COLLAPSE
    << "(org) "
#endif
            << result.first << ", using b * : " << 0.5 << endl;
    outfile << "Collision times: "
#ifdef ORG_COLLAPSE
    << "(org) "
#endif
            << result.first << ", using b * : " << 0.5 << endl;

    int pi_num = 0;
    double error = 0, errorAll = 0, errorALL_org = 0;
    double error_ww = 0, error_org_ww = 0;
#ifdef RS_IN_CELL
    double errorAll_ws = 0, errorALL_org_ws = 0;
#endif

    int skipcounter = 0;
    int all_the_persistent = 0;
    int all_the_persistent_2 = 0;

    for (auto const &it: data.flowLatencies) {
        unsigned int i = it.first;
        auto counterIt = data.flowCounters.find(i);  // 使用 find 而不是 []
        if (counterIt != data.flowCounters.end()) {
            if (counterIt->second > _theta * (data.totalLines / KCU_EPOCH_SIZE + 1)) {
                all_the_persistent++;
                if (counterIt->second > 2 * _theta * (data.totalLines / KCU_EPOCH_SIZE + 1)) {
                    all_the_persistent_2++;
                }
            }
        }
    }

    cout << "\n[DEBUG] Persistence statistics:" << endl;
    cout << "  all_the_persistent (flows above theta threshold): " << all_the_persistent << endl;
    cout << "  all_the_persistent_2 (flows above 2*theta threshold): " << all_the_persistent_2 << endl;
    cout << "  Expected number of persistent flows (theta * unique_flows): "
            << (_theta * data.flowLatencies.size()) << endl;

#ifdef RS_IN_CELL
    double error_array_sum_ws[40] = {};
    double error_array_org_sum_ws[40] = {};
#endif

    double error_array_sum[40] = {};
    double error_array_org_sum[40] = {};

    cout << "\n[DEBUG] Processing bucketElements for error calculation..." << endl;
    for (auto const &element: bucketElements) {
        unsigned int i = element.first;

        auto flowIt = data.flowLatencies.find(i);
        auto counterIt = data.flowCounters.find(i);

        if (flowIt != data.flowLatencies.end() && !flowIt->second.empty() &&
            counterIt != data.flowCounters.end() && counterIt->second > 0) {
            ++pi_num;
#ifdef RS_IN_CELL
            dds_psketch *tmp = dds_psketch_tmp(KCU_type_ptr, i);
            double error_array_ws[40] = {};
            double error_array_org_ws[40] = {};
#endif

            int time = 0;
            double error_array[40] = {};
            double error_array_org[40] = {};

            for (double reqq = 0.6; reqq < 1; reqq += 0.01) {
                double est_ = KCU_QuantileQuery(KCU_type_ptr, i, reqq);

                 double excatVal = query_for_val_vector(flowIt->second, reqq);



                vector<double> sorted_latencies = flowIt->second;
                std::sort(sorted_latencies.begin(), sorted_latencies.end());

                double estQuantile = computePercentile(sorted_latencies, est_);
                double error_org = abs(est_ - excatVal);
                double errorPercentile = abs(estQuantile - reqq);

                error = errorPercentile;

                if (isnan(error)) {
                    error = 1;
                }
                errorAll += error;
                error_array[time] = error;
                error_array_sum[time] += error;
                error_ww += error * counterIt->second / 100; // 使用 counterIt->second

                error = error_org / excatVal;
                if (isnan(error)) {
                    error = 1;
                }
                if (excatVal <= 0) { error = 0; }

                errorALL_org += error;
                error_array_org[time] = error;
                error_array_org_sum[time] += error;
                error_org_ww += error * counterIt->second / 100; // 使用 counterIt->second

#ifdef RS_IN_CELL
                double est = KCU_QuantileQuerywGsamples(tmp, reqq);

                double estQuantile_ws = computePercentile(sorted_latencies, est); // 使用 flowIt->second
                double error_org_ws = abs(est - excatVal);
                double errorPercentile_ws = abs(estQuantile_ws - reqq);

                if (time == 40) {
                    cout << "quantile(with sample) " << reqq << " : " << est_ << endl;
                    cout << "real quantile " << reqq << " : " << excatVal << endl;
                    cout << "est quantile(with sample): " << estQuantile_ws << endl;
                    cout << "Error(with sample): " << errorPercentile_ws << endl;
                }

                error = errorPercentile_ws;
                errorAll_ws += error;
                error_array_ws[time] = error;
                error_array_sum_ws[time] += error;

                error = error_org_ws / excatVal;

                if (excatVal <= 0) { error = 0; }

                errorALL_org_ws += error;
                error_array_org_ws[time] = error;
                error_array_org_sum_ws[time] += error;
#endif
                time++;
            }
#ifdef RS_IN_CELL
            std::cout << "error_array(with sample): ";
            for (int q = 0; q < 40; q++) {
                std::cout << error_array_ws[q] << ",";
                if (q % 10 == 9)std::cout << std::endl;
            }
            std::cout << std::endl;
            std::cout << "error_array_org(with sample): ";
            for (int q = 0; q < 40; q++) {
                std::cout << error_array_org_ws[q] << ",";
                if (q % 10 == 9)std::cout << std::endl;
            }
            std::cout << std::endl;
#endif
        } else {
            skipcounter++;
            bool flowExists = (flowIt != data.flowLatencies.end());
            bool counterExists = (counterIt != data.flowCounters.end());
            bool flowEmpty = flowExists ? flowIt->second.empty() : true;
            bool counterAboveThresh = counterExists ? (counterIt->second > thresh) : false;

            cout << "skip: Flow ID=" << i << " (not monitored: " << flowEmpty
                    << ", exceeded: " << counterAboveThresh << ")" << endl;
        }
    }

    cout << "\n[DEBUG] Final pi_num calculation:" << endl;
    cout << "  pi_num = " << pi_num << endl;
    cout << "  bucketElements size = " << bucketElements.size() << endl;
    cout << "  skipcounter = " << skipcounter << endl;

    if (pi_num > 0) {
        errorAll = errorAll / 40.0 / (double) pi_num;
        error_ww = error_ww / 40.0 / (double) pi_num;
        errorALL_org = errorALL_org / 40.0 / (double) pi_num;
        error_org_ww = error_org_ww / 40.0 / (double) pi_num;
        for (double &array_sum: error_array_sum) {
            array_sum /= pi_num;
        }
        for (double &array_sum_org: error_array_org_sum) {
            array_sum_org /= pi_num;
        }
#ifdef RS_IN_CELL
        errorAll_ws = errorAll_ws / 40.0 / (double) pi_num;
        errorALL_org_ws = errorALL_org_ws / 40.0 / (double) pi_num;
        for (double &array_sum: error_array_sum_ws) {
            array_sum /= pi_num;
        }
        for (double &array_sum_org: error_array_org_sum_ws) {
            array_sum_org /= pi_num;
        }
#endif
    } else {
        cout << "\n[DEBUG] Missing data warning details:" << endl;
        cout << "  bucketElements was " << (bucketElements.empty() ? "EMPTY" : "NOT EMPTY") << endl;
        cout << "  If not empty, but pi_num=0, check conditions:" << endl;
        cout << "    1. data.flowLatencies[i].empty() must be false" << endl;
        cout << "    2. data.flowCounters[i] > 0 must be true" << endl;
        cout << "  Most likely cause: No flow reached max_cell status." << endl;
        cout << "Missing data warning" << endl;
        errorAll = 0;
    }

    cout << "=== Result ===" << endl;
    cout << "all_num: " << all_the_persistent << endl;
    cout << "all_num_2: " << all_the_persistent_2 << endl;
    cout << "est_num: " << pi_num << endl;
    cout << "trace size: " << data.traceSize << endl;
    cout << "Skip: " << skipcounter << endl;
    cout << "RME: " << errorAll << endl;
    cout << "RME_org: " << errorALL_org << endl;
    outfile << "=== Result ===" << endl;
    outfile << "all_num: " << all_the_persistent << endl;
    outfile << "all_num_2: " << all_the_persistent_2 << endl;
    outfile << "est_num: " << pi_num << endl;
    outfile << "trace size: " << data.traceSize << endl;
    outfile << "Skip: " << skipcounter << endl;
    outfile << "RME: " << errorAll << endl;
    outfile << "RME_org: " << errorALL_org << endl;
    std::cout << "error_array_sum: ";
    outfile << "error_array_sum: ";
    for (int q = 0; q < 40; q++) {
        std::cout << error_array_sum[q] << ",";
        outfile << error_array_sum[q] << ",";
    }
    std::cout << std::endl;
    outfile << std::endl;
    std::cout << "error_array_org_sum: ";
    outfile << "error_array_org_sum: ";
    for (int q = 0; q < 40; q++) {
        std::cout << error_array_org_sum[q] << ",";
        outfile << error_array_org_sum[q] << ",";
    }
    std::cout << std::endl;
    outfile << std::endl;

#ifdef RS_IN_CELL
    cout << "=== Result (with sample) ===" << endl;
    cout << "RME (with sample): " << errorAll_ws << endl;
    cout << "RME_org (with sample): " << errorALL_org_ws << endl;
    outfile << "RME (with sample): " << errorAll_ws << endl;
    outfile << "RME_org (with sample): " << errorALL_org_ws << endl;
    std::cout << "error_array_sum (with sample): ";
    outfile << "error_array_sum (with sample): ";
    for (int q = 0; q < 40; q++) {
        std::cout << error_array_sum_ws[q] << ",";
        outfile << error_array_sum_ws[q] << ",";
    }
    std::cout << std::endl;
    outfile << std::endl;
    std::cout << "error_array_org_sum (with sample): ";
    outfile << "error_array_org_sum (with sample): ";
    for (int q = 0; q < 40; q++) {
        std::cout << error_array_org_sum_ws[q] << ",";
        outfile << error_array_org_sum_ws[q] << ",";
    }
    std::cout << std::endl;
    outfile << std::endl;
#endif

    outfile.close();
    string filename_ = "DATA";
    if (filename.find("case") != std::string::npos) {
        filename_ = "CAIDA";
    } else if (filename.find("zipf") != std::string::npos) {
        if (filename.find("_2") != std::string::npos) {
            filename_ = "ZIPF_2";
        } else {
            filename_ = "ZIPF_1";
        }
    }
    if (filename.find("MAWI") != std::string::npos) {
        filename_ = "MAWI";
    }

    int use = KCU_ComputeMem(KCU_type_ptr)[0];
    std::string csv_filename = "../output/quantile/KC_"
                               + filename_ + "mem_" +
                               to_string(use) + "k" +
                               to_string(KCU_type_ptr->k) + ".csv";
    ofstream csv_file(csv_filename);
    if (!csv_file.is_open()) {
        std::cerr << "Error: Cannot open file " << csv_filename << " for writing" << std::endl;
        return -1;
    }

    csv_file << "quantile,percentile_error,relative_error"
#ifdef RS_IN_CELL
            ",percentile_error_ws,relative_error_ws"
#endif
            "\n";

    for (size_t i = 0; i < 40; ++i) {
        csv_file << 0.6 + i * 0.01 << ","
                << error_array_sum[i] << ","
                << error_array_org_sum[i]
#ifdef RS_IN_CELL
                << ","
                   << error_array_sum_ws[i] << ","
                   << error_array_org_sum_ws[i]
#endif
                << "\n";
    }
    csv_file.close();

    csv_filename = "../output/sum_" + filename_ + ".csv";
    ofstream csv_file_2(csv_filename, ios_base::app);
    if (!csv_file_2.is_open()) {
        std::cerr << "Error: Cannot open file " << csv_filename << " for writing" << std::endl;
        return -1;
    }

    csv_file_2 << use / pow(10.0, 6.0) << "," << pi_num << "," << errorAll << "," << errorALL_org
#ifdef RS_IN_CELL
    << "," << errorAll_ws << ","
       << errorALL_org_ws
#endif
            << "\n";

    csv_file_2.close();
    std::cout << "[INFO] CSV data saved to " << csv_filename << std::endl;

    csv_filename = "../output/sum_ww_" + filename_ + ".csv";
    ofstream csv_file_3(csv_filename, ios_base::app);
    if (!csv_file_3.is_open()) {
        std::cerr << "Error: Cannot open file " << csv_filename << " for writing" << std::endl;
        return -1;
    }

    csv_file_3 << use / pow(10.0, 6.0) << "," << error_ww << "," << error_org_ww
#ifdef ORG_COLLAPSE
    << "," << 2
#else
            << "," << 1
#endif
            << "\n";

    csv_file_3.close();
    std::cout << "[INFO] CSV data saved to " << csv_filename << std::endl;

    return 0;
}

// Original Tail_letency_quantile_test function as a wrapper
int Tail_letency_quantile_test(std::string filename, KCU_type *KCU_type_ptr) {
    TraceData data = readTraceFile(filename);
    if (data.totalLines == 0) return -1;

    populateKCUSketch(KCU_type_ptr, data);

    return Tail_latency_quantile_test_with_data(filename, KCU_type_ptr, data);
}

// New function for running multiple tests with same data
int run_multiple_tests_with_same_data(const std::string &filename,
                                      const vector<KCU_type *> &kcu_sketches) {
    TraceData data = readTraceFile(filename);
    if (data.totalLines == 0) return -1;

    for (size_t i = 0; i < kcu_sketches.size(); i++) {
        cout << "\n=== Running test " << (i + 1) << " of " << kcu_sketches.size() << " ===" << endl;

        TraceData tempData = data;

        populateKCUSketch(kcu_sketches[i], tempData);

        Tail_latency_quantile_test_with_data(filename, kcu_sketches[i], tempData);
    }

    return 0;
}

void processCounters(
    std::unordered_map<uint32_t, uint32_t> &countersMap,
    std::map<uint32_t, uint32_t> &bucketElements, double &AAE, string filename, KCU_type *kcu) {
    std::unordered_map<uint32_t, uint32_t> minheap1;
    std::unordered_map<uint32_t, uint32_t> minheap2;
    std::unordered_map<uint32_t, uint32_t> minheap3;
    // sort
    std::vector<std::pair<uint32_t, uint32_t> > elements;
    elements.reserve(countersMap.size());

    for (auto &pair: countersMap) {
        elements.emplace_back(pair);
    }


    std::sort(elements.begin(), elements.end(),
              [](const std::pair<uint32_t, uint32_t> &a,
                 const std::pair<uint32_t, uint32_t> &b) {
                  return a.second > b.second;
              });


    minheap1.clear();
    size_t count1 = std::min(elements.size(), bucketElements.size());
    for (size_t i = 0; i < count1; ++i) {
        minheap1[elements[i].first] = elements[i].second;
    }

    minheap2.clear();
    size_t count2 = std::min(elements.size(), bucketElements.size() * 2);
    for (size_t i = 0; i < count2; ++i) {
        minheap2[elements[i].first] = elements[i].second;
    }
    minheap3.clear();
    size_t count3 = std::min(elements.size(), bucketElements.size() * 4);
    for (size_t i = 0; i < count3; ++i) {
        minheap3[elements[i].first] = elements[i].second;
    }

    // count
    int countInHeap1 = 0, countInHeap2 = 0, countInHeap3 = 0;

    for (const auto &pair: bucketElements) {
        uint32_t key = pair.first;


        if (minheap1.find(key) != minheap1.end()) {
            countInHeap1++;
        }
        if (minheap2.find(key) != minheap1.end()) {
            countInHeap2++;
        }
        if (minheap3.find(key) != minheap1.end()) {
            countInHeap3++;
        }
        if (countersMap.find(key) != countersMap.end()) {
            AAE += abs((int) countersMap[key] - (int) pair.second);
        }
    }

    AAE /= (double) bucketElements.size();
    std::cout << "Result：" << std::endl;
    std::cout << "key size: " << elements.size() << std::endl;
    std::cout << "key in " << bucketElements.size() << "-top persistence: " << countInHeap1 << "Rate: " << (double)
            countInHeap1 / (double) bucketElements.size() << std::endl;
    std::cout << "key in " << bucketElements.size() * 2 << "-top persistence: " << countInHeap2 << "Rate: " << (double)
            countInHeap2 / (double) bucketElements.size() / 2.0 << std::endl;
    std::cout << "AAE: " << AAE << std::endl;


    string filename_ = "DATA";;
    if (filename.find("case") != std::string::npos) {
        filename_ = "CAIDA";
    } else if (filename.find("zipf") != std::string::npos) {
        if (filename.find("_2") != std::string::npos) {
            filename_ = "ZIPF_2";
        } else {
            filename_ = "ZIPF_1";
        }
    }

    if (filename.find("MAWI") != std::string::npos) {
        filename_ = "MAWI";
    }
    std::string filename_csv = "../" + static_cast<std::string>(filename_) +
                               "Persistence_QUARTZ.csv";

    std::ofstream outfile(filename_csv, ios::app);

    if (!outfile.is_open()) {
        std::cerr << "Error: Cannot open file " << filename_csv << " for writing" << std::endl;
        return;
    }

    outfile << bucketElements.size() << "," << AAE << "," << (double)
            countInHeap1 / (double) bucketElements.size() << "\n";

    outfile.close();
    std::cout << "[INFO] QUARTZ persistence results saved to " << filename << std::endl;
}


int Persistence_filter_Test(std::string filename) {
    int nonZerosLatency = 0;
    int traceSize = 0;
    unordered_map<uint32_t, vector<double> > map; // key (flow id), vector of latencies of the flow id
    unordered_map<uint32_t, uint32_t> countersMap; // key (flow id) , counter
    unordered_map<uint32_t, uint32_t> lastUpdatedEpoch;
    vector<pair<uint32_t, double> > vec;
    cout << "debug: " << filename << endl;

    //Parse the data
    ifstream infile;
    string line;
    uint32_t ID;
    double latency;
    infile.open(filename);

    if (!infile.is_open()) {
        cout << "can't open file " << filename << endl;
        return -1;
    }

    int lineCount = 0;
    while (!infile.eof()) {
        std::vector<std::string> splits;
        std::string split;

        // Read a line
        getline(infile, line);
        if (line.empty()) continue;

        lineCount++;
        if (lineCount % 100000 == 0) {
            cout << "processed " << lineCount << " lines" << endl;
        }

        std::istringstream iss(line);

        char delim = ',';
        if (filename.find("case") != std::string::npos)delim = ' ';
        while (std::getline(iss, split, delim)) {
            splits.push_back(split);
        }

        try {
            ID = (uint32_t) stoull(splits[0]);
            latency = stoi(splits[1]);
        } catch (const std::exception &e) {
            cout << "line: " << lineCount << " transform error " << e.what() << ", content: " << line << endl;
            continue;
        }

        if (latency > 0) {
            ++nonZerosLatency;
        } else {
            continue;
        }

        int epoch = lineCount / KCU_EPOCH_SIZE + 1;

        map[ID].push_back(latency);
        auto it = lastUpdatedEpoch.find(ID);
        if (it == lastUpdatedEpoch.end() || it->second != epoch) {
            ++countersMap[ID];
            lastUpdatedEpoch[ID] = epoch;
        }
        ++traceSize;
        vec.emplace_back(ID, latency);
    }

    infile.close();


    for (int i = 500; i <= 1500; i += 100) {

#ifdef USE_IMPROVED_VERSION
        i /= SKETCH_IN_BUCKET;
#endif

        KCU_type *KCU = KCU_Init(1.0 / i, 0.05, 0, i * 100);
        for (pair<uint32_t, double> p: vec) {
            KCU_UpdateLatency(KCU, p.first, p.second);
        }
        std::map<uint32_t, uint32_t> bucketElements = KCU_Output(KCU);
        double AAE_P = 0;
        processCounters(countersMap, bucketElements, AAE_P, filename, KCU);
        std::map<uint32_t, uint32_t> bucketElements1 = KCU_Output3(KCU);
        AAE_P = 0;
        processCounters(countersMap, bucketElements1, AAE_P, filename, KCU);
    }

    return 0;
}


//####################################################################

int Memory_THP_test(std::string filename, KCU_type *sketch, double alpha) {
    int traceSize = 0;
    unordered_map<uint32_t, vector<double> > map; // key (flow id), vector of latencies of the flow id
    unordered_map<uint32_t, uint32_t> countersMap; // key (flow id) , counter

    //Parse the data
    ifstream infile;
    string line;
    uint32_t ID;
    double latency;
    infile.open(filename);

    while (!infile.eof()) {
        std::vector<std::string> splits;
        std::string split;

        // Read a line
        getline(infile, line);
        if (line.empty()) continue;
        std::istringstream iss(line);

        char delim = ',';
        if (filename.find("case") != std::string::npos)delim = ' ';
        while (std::getline(iss, split, delim)) {
            splits.push_back(split);
        }

        ID = (uint32_t) stoull(splits[0]);
        latency = stoi(splits[1]);

        map[ID].push_back(latency);
        ++countersMap[ID];
        ++traceSize;
    }
    double left = 0.18,right = 0.5,step = 10;
    //double left = 0.4,right = 1.1,step = 10;
  //  double left = 0.9, right = 1.8, step = 10;
    int times = 0;
    int k_ = sketch->k;
    for (double l = left; l < right; l += (right - left) / step) {
        times++;
        int max_buckets_for_per = CalculateMaxBucketsForMemory((int) (l * 1000000), k_);
        // if (times > 7)alpha = 0.005;
        KCU_type *KCU = KCU_Init((double) 1 / (double) k_, alpha, 0, max_buckets_for_per);

        clock_t begint, endt;
        struct timeb begintb, endtb;
        double time;
        unsigned int numItems = 0;

        begint = clock();
        ftime(&begintb);
        for (auto pair: map) {
            for (auto ltc: pair.second) {
                KCU_UpdateLatency(KCU, pair.first, ltc);
                ++numItems;
            }
        }

        endt = clock();
        ftime(&endtb);
        time = ((double) (endt - begint)) / CLK_PER_SEC;


        string filename_ = "DATA";;
        if (filename.find("case") != std::string::npos) {
            filename_ = "CAIDA";
        } else if (filename.find("zipf") != std::string::npos) {
            if (filename.find("_2") != std::string::npos) {
                filename_ = "ZIPF_2";
            } else {
                filename_ = "ZIPF_1";
            }
        }
        if (filename.find("MAWI") != std::string::npos) {
            filename_ = "MAWI";
        }
        std::vector<int> mem = KCU_ComputeMem(KCU);
        std::string csv_filename = "../KC_THP_" + filename_ + "b_" + to_string(KCU->k) + ".csv";
        // output to csv
        std::ofstream csv_file(csv_filename, ios::app);
        if (!csv_file.is_open()) {
            std::cerr << "Error: Cannot open file " << csv_filename << " for writing" << std::endl;
            return -1;
        }
        double size = 632489;
        if (filename.find("MAWI") != string::npos) {
            size = 588529;
        }
        if (filename.find("_1") != string::npos) {
            size = 338859;
        }

        if (filename.find("_2") != string::npos) {
            size = 342104;
        }
        csv_file << size / time << "," << mem[0] / pow(10.0, 6.0) << "," << mem[1]
                << "," << mem[2] << "," << mem[3] << "," << mem[4] << "\n";
        csv_file.close();
        std::cout << "[INFO] CSV data saved to " << csv_filename << std::endl;

        printf("Update_KC_UDD %d pairs took %lfs %f alpha\n", numItems, time, alpha);
        KCU_Destroy(KCU);
    }
    return 0;
}

//##################################################################

int main(int argc, char *argv[]) {
    if (argc != 6) {
        std::cout << "arguments expected:" << std::endl;
        std::cout << "(1) file name" << std::endl;
        std::cout << "(2) theta " << std::endl;
        std::cout << "(3) alpha" << std::endl;
        std::cout << "(4) max buckets" << std::endl;
        std::cout << "(5) type" << std::endl;
        return -1;
    }


    _theta = atof(argv[2]);
    alpha = atof(argv[3]);
    max_buckets = atoi(argv[4]);
    int type = atoi(argv[5]);

    //theta = (double) theta / (double) (1+sqrt(epsilon));

    long long bucket_size = (int) ((double) (pow(_theta, -1)));
    long long samples = (int) (pow(alpha, -1) / (double) (1 + KCU_CELL_SIZE));


    KCU_type *KCU = KCU_Init((double) 1 / (double) bucket_size, alpha, samples, max_buckets);
    int result = -1;
    //double left = 0.18,right = 0.5,step = 10;
   // double left = 0.4,right = 1.1,step = 10;
    double left = 0.9, right = 1.8, step = 10;
 // double left = 0.18,right = 0.5,step = 10;

    switch (type) {
        case 1:
            result = Memory_THP_test(argv[1], KCU, alpha);
            break;
        case 2:
            result = Persistence_filter_Test(argv[1]);
            break;
        case 3:
            result = Tail_letency_quantile_test(argv[1], KCU);
            break;
        case 4:
            int times = 0;
            TraceData data = readTraceFile(argv[1]);
            if (data.totalLines == 0) return -1;


            for (double hajime = left; hajime < right; hajime += (right - left) / step) {
                times++;
                int max_buckets_for_per = CalculateMaxBucketsForMemory((int) (hajime * 1000000), bucket_size);
                if (times > 7)alpha = 0.005;
                KCU_type *mine = KCU_Init((double) 1 / (double) bucket_size, alpha, samples, max_buckets_for_per);
                populateKCUSketch(mine, data);

                 result = Tail_latency_quantile_test_with_data(argv[1], mine, data);

                KCU_Destroy(mine);
            }
            break;
    }


    return result;
}
