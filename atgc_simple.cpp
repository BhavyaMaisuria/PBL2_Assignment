#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <string>
#include <vector>

// UI Helpers
const std::string RESET = "\033[0m";
const std::string BOLD = "\033[1m";
const std::string CYAN = "\033[96m";
const std::string GREEN = "\033[92m";
const std::string YELLOW = "\033[93m";
const std::string RED = "\033[91m";
const std::string PURPLE = "\033[95m";

void printHeader(const std::string& title) {
    std::cout << "\n" << CYAN << BOLD 
              << "========================================\n"
              << "  " << title << "\n"
              << "========================================" 
              << RESET << "\n";
}

void loadGenomicData(std::vector<std::string> &species,
                     std::vector<std::string> &fragments) {
    std::ifstream file("dna_data.txt");
    if (!file) {
        std::cout << RED << "❌ Error: dna_data.txt not found in directory!" << RESET << "\n";
        return;
    }

    std::string line;
    bool isFragment = false;
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty() || line[0] == '>')
            continue; // Skip empty lines and FASTA headers

        if (line[0] == '#') {
            if (line.find("FRAGMENT") != std::string::npos) {
                isFragment = true;
            }
            continue;
        }

        if (!isFragment) {
            species.push_back(line);
        } else {
            fragments.push_back(line);
        }
    }
    std::cout << GREEN << "✅ Database Loaded: " << BOLD << species.size() << RESET << GREEN << " DNA Strands, " 
              << BOLD << fragments.size() << RESET << GREEN << " Fragments. 🧬\n" << RESET;
}

void runAlignment(const std::string& s1, const std::string& s2, bool isGlobal) {
    printHeader("🔍 Alignment Analysis");
    std::cout << PURPLE << "🧬 Seq A: " << RESET << s1 << "\n"
              << PURPLE << "🎯 Seq B: " << RESET << s2 << "\n\n";

    int m = s1.length(), n = s2.length();
    std::vector<std::vector<int>> dp(m + 1, std::vector<int>(n + 1, 0));
    int match = 2, mismatch = -1, gap = -2;

    if (isGlobal) {
        for (int i = 0; i <= m; i++) dp[i][0] = i * gap;
        for (int j = 0; j <= n; j++) dp[0][j] = j * gap;
    }

    int maxS = 0;
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int score = (s1[i - 1] == s2[j - 1]) ? match : mismatch;
            int res = std::max({dp[i - 1][j - 1] + score,
                                dp[i - 1][j] + gap,
                                dp[i][j - 1] + gap});
            dp[i][j] = isGlobal ? res : std::max(0, res);
            if (dp[i][j] > maxS) maxS = dp[i][j];
        }
    }

    int finalScore = isGlobal ? dp[m][n] : maxS;
    float sim = (std::max(m, n) > 0)
                    ? (float)finalScore / (std::max(m, n) * match) * 100
                    : 0;
    
    if (sim < 0) sim = 0.0f;
    if (sim > 100) sim = 100.0f;

    std::cout << CYAN << "▶ Mode:       " << RESET << (isGlobal ? "🌍 Global" : "📍 Local") << "\n"
              << CYAN << "▶ Score:      " << RESET << BOLD << finalScore << RESET << "\n"
              << CYAN << "▶ Similarity: " << RESET << BOLD << std::fixed << std::setprecision(2) << sim << "%" << RESET << "\n";

    std::cout << CYAN << "▶ Status:     " << RESET;
    if (sim > 90)      std::cout << GREEN << BOLD << "🟢 CONSERVED SEQUENCE\n" << RESET;
    else if (sim > 50) std::cout << YELLOW << BOLD << "🟡 DIVERGENT EVOLUTION\n" << RESET;
    else               std::cout << RED << BOLD << "🔴 HIGH MUTATION\n" << RESET;
    std::cout << CYAN << "----------------------------------------\n" << RESET;
}

void runFragmentAssembly(const std::vector<std::string>& frags) {
    printHeader("🧩 Fragment Assembly");
    for (size_t i = 0; i < frags.size(); i++)
        std::cout << PURPLE << "🔬 Node[" << i << "]: " << RESET << frags[i] << "\n";
    std::cout << "\n";

    int n = frags.size();
    std::vector<int> inDegree(n, 0);
    std::vector<std::vector<int>> adj(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j && frags[i].length() >= 2 && frags[j].length() >= 2) {
                if (frags[i].substr(frags[i].length() - 2) == frags[j].substr(0, 2)) {
                    adj[i].push_back(j);
                    inDegree[j]++;
                }
            }
        }
    }

    std::queue<int> q;
    for (int i = 0; i < n; i++)
        if (inDegree[i] == 0) q.push(i);

    std::string finalDNA = "";
    int processedCount = 0;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        processedCount++;
        finalDNA += (finalDNA == "" ? frags[u] : frags[u].substr(2));
        for (int v : adj[u])
            if (--inDegree[v] == 0) q.push(v);
    }

    if (processedCount < n) {
        std::cout << RED << "⚠️ Warning: Fragment graph contains cycles or disjoint sets.\n" << RESET;
    }

    std::cout << CYAN << "▶ Target Assembly: " << GREEN << BOLD 
              << (finalDNA.empty() ? "N/A" : finalDNA) << RESET << "\n"
              << CYAN << "----------------------------------------\n" << RESET;
}

int main() {
    std::cout << "\n" << CYAN << BOLD 
              << "========================================\n"
              << "       🌟 GENOMATCH PRO BOOTING... 🌟\n"
              << "========================================" << RESET << "\n";
    
    std::vector<std::string> species, fragments;
    loadGenomicData(species, fragments);

    int choice;
    while (true) {
        std::cout << "\n" << YELLOW << BOLD << "🚀 MAIN MENU\n" << RESET
                  << "  🌍  1. Global Alignment\n"
                  << "  📍  2. Local Alignment\n"
                  << "  ⚙️  3. Fragment Assembly\n"
                  << "  ❌  0. Exit\n"
                  << YELLOW << "\nCommand > " << RESET;
        
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::cin.ignore(10000, '\n');
            continue;
        }

        if (choice == 0) {
            std::cout << CYAN << "Goodbye! 👋\n" << RESET;
            break;
        }

        if (choice == 1) {
            if (species.size() >= 2) runAlignment(species[0], species[1], true);
            else std::cout << RED << "⚠️ Not enough DNA Strands for Global Alignment.\n" << RESET;
        } else if (choice == 2) {
            if (species.size() >= 4) runAlignment(species[2], species[3], false);
            else if (species.size() >= 2) runAlignment(species[0], species[1], false);
            else std::cout << RED << "⚠️ Not enough DNA Strands for Local Alignment.\n" << RESET;
        } else if (choice == 3) {
            if (!fragments.empty()) runFragmentAssembly(fragments);
            else std::cout << RED << "⚠️ No Fragments available for assembly.\n" << RESET;
        } else {
            std::cout << RED << "⚠️ Invalid Command. Try again!\n" << RESET;
        }
    }
    return 0;
}
