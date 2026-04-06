#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <queue>
#include <thread>
#include <chrono>
#include <atomic>
#include <iomanip>

using namespace std;

// ==========================================
// 1. UI THEMES & MULTI-THREADING GLOBALS
// ==========================================
struct Theme { string p, s, a, name; };
Theme bio = {"\033[96m", "\033[95m", "\033[92m", "BIO-CYAN 🧬"};
Theme lava = {"\033[91m", "\033[93m", "\033[94m", "LAVA-RED 🔥"};
Theme cur = bio;

atomic<int> pulse(0);
atomic<bool> active(true);
const vector<string> COLORS = {"\033[94m", "\033[96m", "\033[92m", "\033[95m", "\033[93m"};

void pulseWorker() {
    while (active) {
        this_thread::sleep_for(chrono::seconds(2));
        pulse = (pulse + 1) % COLORS.size();
    }
}

// ==========================================
// 2. DATA LOADER (FILE I/O)
// ==========================================
void loadGenomicData(vector<string>& species, vector<string>& fragments) {
    ifstream file("dna_data.txt");
    if (!file) { 
        cout << "\033[31m ❌ Error: dna_data.txt not found in directory!\033[0m" << endl; 
        return; 
    }
    
    string line;
    while (getline(file, line)) {
        if (line.empty() || line[0] == '#') continue; 
        if (line[0] == '>') continue; // Skip FASTA headers
        
        if (line.length() > 10) species.push_back(line);
        else fragments.push_back(line);
    }
    cout << "\033[92m ✅ Database Sync Complete: " << species.size() << " DNA Strands, " << fragments.size() << " Fragments.\033[0m\n" << endl;
}

// ==========================================
// 3. ALGORITHM: DP ALIGNMENT (GLOBAL & LOCAL)
// ==========================================
void runAlignment(string s1, string s2, bool isGlobal) {
    // A. Echo Input
    cout << cur.p << "\n[📡 LIVE DATA FEED]" << endl;
    cout << " » Seq A: \033[0m" << s1 << endl;
    cout << cur.p << " » Seq B: \033[0m" << s2 << endl;
    cout << cur.s << "------------------------------------------------" << "\033[0m" << endl;

    // B. Algorithm Core (Needleman-Wunsch / Smith-Waterman)
    int m = s1.length(), n = s2.length();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    int match = 2, mismatch = -1, gap = -2;

    if (isGlobal) { 
        for (int i = 0; i <= m; i++) dp[i][0] = i * gap;
        for (int j = 0; j <= n; j++) dp[0][j] = j * gap;
    }

    int maxS = 0;
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int score = (s1[i - 1] == s2[j - 1]) ? match : mismatch;
            int res = max({dp[i - 1][j - 1] + score, dp[i - 1][j] + gap, dp[i][j - 1] + gap});
            dp[i][j] = isGlobal ? res : max(0, res);
            if (dp[i][j] > maxS) maxS = dp[i][j];
        }
    }

    // C. Decode Multi-Info Output
    int finalScore = isGlobal ? dp[m][n] : maxS;
    float sim = (float)finalScore / (max(m, n) * match) * 100;
    if (sim < 0) sim = 0.0;
    if (sim > 100) sim = 100.0;

    cout << cur.s << "[📊 DECODED ANALYSIS REPORT]" << endl;
    cout << "\033[0m » Target Paradigm: " << (isGlobal ? "Global (Needleman-Wunsch / Evolutionary)" : "Local (Smith-Waterman / Mutation)") << endl;
    cout << " » Match Score:     " << cur.a << finalScore << "\033[0m" << endl;
    cout << " » Similarity:      " << fixed << setprecision(2) << sim << "%" << endl;
    
    if (sim > 90) cout << "\033[92m » System Status:   CONSERVED SEQUENCE ✅\033[0m" << endl;
    else if (sim > 50) cout << "\033[93m » System Status:   DIVERGENT EVOLUTION 🧬\033[0m" << endl;
    else cout << "\033[91m » System Status:   HIGH MUTATION ALERT ⚠️\033[0m" << endl;
    cout << cur.s << "------------------------------------------------" << "\033[0m" << endl;
}

// ==========================================
// 4. ALGORITHM: GRAPH ASSEMBLY (KAHN'S DAG)
// ==========================================
void runFragmentAssembly(vector<string> frags) {
    // A. Echo Input
    cout << cur.p << "\n[📡 FRAGMENT FEED]" << "\033[0m" << endl;
    for(size_t i=0; i<frags.size(); i++) cout << " » Node[" << i << "]: " << frags[i] << endl;
    cout << cur.s << "------------------------------------------------" << "\033[0m" << endl;

    // B. Build Graph & In-Degree Array
    int n = frags.size();
    vector<int> inDegree(n, 0);
    vector<vector<int>> adj(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j && frags[i].substr(frags[i].length()-2) == frags[j].substr(0,2)) {
                adj[i].push_back(j);
                inDegree[j]++;
            }
        }
    }

    // C. Kahn's Topological Sort
    queue<int> q;
    for (int i = 0; i < n; i++) if (inDegree[i] == 0) q.push(i);

    cout << cur.s << "[🧩 RECONSTRUCTION PROGRESS]" << "\033[0m" << endl;
    string finalDNA = "";
    while (!q.empty()) {
        int u = q.front(); q.pop();
        cout << " » Processing: " << frags[u] << " ➔ ";
        finalDNA += (finalDNA == "" ? frags[u] : frags[u].substr(2));
        for (int v : adj[u]) if (--inDegree[v] == 0) q.push(v);
    }
    cout << "DONE" << endl;
    cout << " » Target Assembly: \033[92m" << finalDNA << "\033[0m" << endl;
    cout << cur.s << "------------------------------------------------" << "\033[0m" << endl;
}

// ==========================================
// 5. MAIN MENU & EXECUTION LOOP
// ==========================================
void renderMenu() {
    string p = COLORS[pulse];
    #ifdef _WIN32
        system("cls");
    #else
        system("clear");
    #endif
    cout << p << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━" << endl;
    cout << "   🧬 GENOMATCH PRO | ACTIVE THEME: " << cur.name << endl;
    cout << "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\033[0m" << endl;
    cout << cur.p << "  1. 🌐 Global Alignment (Evolution Tracker)" << endl;
    cout << "  2. 📍 Local Alignment (Mutation Scanner)" << endl;
    cout << "  3. 🔗 Fragment Assembly (Kahn's DAG Algorithm)" << endl;
    cout << "  4. 🎨 Switch UI Theme" << endl;
    cout << "  0. ❌ Shutdown System" << endl;
    cout << "\n " << p << " Command > \033[0m";
}

int main() {
    // Start Background Thread for UI
    thread t(pulseWorker);
    
    // Auto-Load Data on Startup
    vector<string> species, fragments;
    cout << "\n\033[96m[SYSTEM] Booting GenoMatch Engine...\033[0m" << endl;
    this_thread::sleep_for(chrono::seconds(1)); 
    loadGenomicData(species, fragments);
    this_thread::sleep_for(chrono::seconds(1));

    int choice;
    while (true) {
        renderMenu();
        if (!(cin >> choice)) { // Handle invalid string inputs
            cin.clear(); cin.ignore(1000, '\n'); continue;
        }

        if (choice == 0) break;
        
        // Execute based on choice
        if (choice == 1 && species.size() >= 2) {
            runAlignment(species[0], species[1], true);
        } 
        else if (choice == 2 && species.size() >= 2) {
            // Using a subset of the string to simulate a local signature search
            runAlignment(species[0], species[2], false); 
        } 
        else if (choice == 3 && fragments.size() >= 3) {
            runFragmentAssembly(fragments);
        } 
        else if (choice == 4) {
            cur = (cur.name == bio.name) ? lava : bio;
            continue; // Skip the pause
        } else {
            cout << "\033[31m ⚠️  Invalid Command or Data Missing.\033[0m" << endl;
        }

        cout << "\n\033[90mPress Enter to sync with Command Center...\033[0m";
        cin.ignore(); cin.get();
    }

    // Clean Exit
    active = false;
    t.join();
    return 0;
}