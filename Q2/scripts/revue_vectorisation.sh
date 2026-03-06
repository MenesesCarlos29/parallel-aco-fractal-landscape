#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./scripts/revue_vectorisation.sh <fichier_cpp>

Description:
  Compile un seul fichier C++ de Q2 avec le rapport de vectorisation GCC
  et genere un resume sous forme de table dans Q2/results/reviews/.

Exemples:
  ./scripts/revue_vectorisation.sh ant.cpp
  ./scripts/revue_vectorisation.sh ant_simu.cpp
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if [[ $# -ne 1 ]]; then
    usage
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
Q2_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
INPUT_PATH="$1"

if [[ "${INPUT_PATH}" == Q2/* ]]; then
    INPUT_PATH="${INPUT_PATH#Q2/}"
fi

if [[ "${INPUT_PATH}" == /* ]]; then
    TARGET_PATH="${INPUT_PATH}"
else
    TARGET_PATH="${Q2_DIR}/${INPUT_PATH}"
fi

if [[ "${TARGET_PATH}" != *.cpp && -f "${TARGET_PATH}.cpp" ]]; then
    TARGET_PATH="${TARGET_PATH}.cpp"
fi

TARGET_REAL="$(realpath "${TARGET_PATH}" 2>/dev/null || true)"
if [[ -z "${TARGET_REAL}" || ! -f "${TARGET_REAL}" ]]; then
    echo "Erreur: fichier introuvable: ${INPUT_PATH}"
    exit 2
fi

if [[ "${TARGET_REAL}" != "${Q2_DIR}/"* ]]; then
    echo "Erreur: ce script accepte uniquement des fichiers dans Q2/."
    exit 3
fi

if [[ "${TARGET_REAL}" != *.cpp ]]; then
    echo "Erreur: fichier non supporte. Donne un .cpp a analyser."
    exit 4
fi

TARGET_REL="${TARGET_REAL#${Q2_DIR}/}"
TARGET_BASE="$(basename "${TARGET_REL}")"
TARGET_STEM="${TARGET_BASE%.cpp}"

OUT_DIR="${Q2_DIR}/results/reviews"
mkdir -p "${OUT_DIR}"

TMP_DIR="$(mktemp -d)"
RAW_REPORT="${TMP_DIR}/vectorisation_brut.txt"
CSV_REPORT="${OUT_DIR}/revue_vectorisation_${TARGET_STEM}.csv"
CSV_FOCUS="${OUT_DIR}/revue_vectorisation_${TARGET_STEM}_focus.csv"
MD_REPORT="${OUT_DIR}/revue_vectorisation_${TARGET_STEM}.md"
OBJ_FILE="${TMP_DIR}/${TARGET_STEM}.o"

CXX_BIN="${CXX:-g++}"

pushd "${Q2_DIR}" >/dev/null
if ! env -u DEBUG "${CXX_BIN}" \
    -std=c++17 -O3 -march=native -Wall \
    -fopt-info-vec-optimized \
    -fopt-info-vec-missed \
    -c "${TARGET_REL}" \
    -o "${OBJ_FILE}" 2>"${RAW_REPORT}"; then
    popd >/dev/null
    echo "Erreur: echec de compilation pour ${TARGET_REL}."
    echo "Sortie compilateur:"
    cat "${RAW_REPORT}"
    rm -rf "${TMP_DIR}"
    exit 5
fi
popd >/dev/null

awk -v rel="${TARGET_REL}" -v base="${TARGET_BASE}" '
BEGIN {
    print "ligne,colonne,statut,message";
}
{
    if (match($0, /^([^:]+):([0-9]+):([0-9]+): (optimized|missed): (.*)$/, m)) {
        src = m[1];
        line = m[2];
        col = m[3];
        status = m[4];
        msg = m[5];
        if (src == rel || src == base) {
            gsub(/"/, "'\''", msg);
            print line "," col "," status ",\"" msg "\"";
        }
    }
}
' "${RAW_REPORT}" > "${CSV_REPORT}"

awk -F, '
NR == 1 {
    print $0;
    next;
}
{
    status = $3;
    msg = $4;
    sub(/^"/, "", msg);
    sub(/"$/, "", msg);

    if (msg ~ /statement can throw an exception/) next;
    if (msg ~ /resx [0-9]+/) next;

    if (msg ~ /vectorize loop/ ||
        msg ~ /not vectorized: loop/ ||
        msg ~ /not vectorized: number of iterations/ ||
        msg ~ /not vectorized: unsupported use/ ||
        msg ~ /not vectorized: control flow/ ||
        msg ~ /not vectorized: complicated access pattern/) {
        print $0;
    }
}
' "${CSV_REPORT}" > "${CSV_FOCUS}"

TOTAL_SIGNAUX="$(awk 'NR > 1 { c++ } END { print c + 0 }' "${CSV_REPORT}")"
NB_OPT="$(awk -F, 'NR > 1 && $3 == "optimized" { c++ } END { print c + 0 }' "${CSV_REPORT}")"
NB_MISSED="$(awk -F, 'NR > 1 && $3 == "missed" { c++ } END { print c + 0 }' "${CSV_REPORT}")"
TOTAL_FOCUS="$(awk 'NR > 1 { c++ } END { print c + 0 }' "${CSV_FOCUS}")"
FOCUS_OPT="$(awk -F, 'NR > 1 && $3 == "optimized" { c++ } END { print c + 0 }' "${CSV_FOCUS}")"
FOCUS_MISSED="$(awk -F, 'NR > 1 && $3 == "missed" { c++ } END { print c + 0 }' "${CSV_FOCUS}")"

TOP_MISSED="${TMP_DIR}/top_missed.txt"
awk -F, '
NR > 1 && $3 == "missed" {
    msg = $4;
    sub(/^"/, "", msg);
    sub(/"$/, "", msg);
    count[msg]++;
}
END {
    for (k in count) {
        print count[k] "\t" k;
    }
}
' "${CSV_REPORT}" | sort -nr > "${TOP_MISSED}"

TOP_MISSED_FOCUS="${TMP_DIR}/top_missed_focus.txt"
awk -F, '
NR > 1 && $3 == "missed" {
    msg = $4;
    sub(/^"/, "", msg);
    sub(/"$/, "", msg);
    count[msg]++;
}
END {
    for (k in count) {
        print count[k] "\t" k;
    }
}
' "${CSV_FOCUS}" | sort -nr > "${TOP_MISSED_FOCUS}"

{
    echo "# Revue vectorisation Q2"
    echo
    echo "- Fichier analyse: \`${TARGET_REL}\`"
    echo "- Signaux detectes: **${TOTAL_SIGNAUX}**"
    echo "- Signaux \`optimized\`: **${NB_OPT}**"
    echo "- Signaux \`missed\`: **${NB_MISSED}**"
    echo "- Signaux utiles (loops/causes actionnables): **${TOTAL_FOCUS}**"
    echo "- Utiles \`optimized\`: **${FOCUS_OPT}**"
    echo "- Utiles \`missed\`: **${FOCUS_MISSED}**"
    echo
    echo "## Table prioritaire (a corriger en premier)"
    echo
    echo "| Ligne | Colonne | Statut | Message |"
    echo "|---:|---:|:---:|---|"

    if [[ "${TOTAL_FOCUS}" -eq 0 ]]; then
        echo "| - | - | - | Aucun signal prioritaire detecte |"
    else
        awk -F, '
        NR > 1 {
            line = $1;
            col = $2;
            status = $3;
            msg = $4;
            sub(/^"/, "", msg);
            sub(/"$/, "", msg);
            gsub(/\|/, "\\|", msg);
            print "| " line " | " col " | " status " | " msg " |";
        }
        ' "${CSV_FOCUS}"
    fi

    echo
    echo "## Causes principales (prioritaires)"
    echo
    echo "| Occurrences | Cause |"
    echo "|---:|---|"

    if [[ "${FOCUS_MISSED}" -eq 0 ]]; then
        echo "| 0 | Aucune cause missed prioritaire |"
    else
        awk -F'\t' '
        NR <= 10 {
            gsub(/\|/, "\\|", $2);
            print "| " $1 " | " $2 " |";
        }
        ' "${TOP_MISSED_FOCUS}"
    fi

    echo
    echo "## Table complete (tous signaux)"
    echo
    echo "| Ligne | Colonne | Statut | Message |"
    echo "|---:|---:|:---:|---|"

    if [[ "${TOTAL_SIGNAUX}" -eq 0 ]]; then
        echo "| - | - | - | Aucun signal pour ce fichier |"
    else
        awk -F, '
        NR > 1 {
            line = $1;
            col = $2;
            status = $3;
            msg = $4;
            sub(/^"/, "", msg);
            sub(/"$/, "", msg);
            gsub(/\|/, "\\|", msg);
            print "| " line " | " col " | " status " | " msg " |";
        }
        ' "${CSV_REPORT}"
    fi

    echo
    echo "## Causes principales (complet)"
    echo
    echo "| Occurrences | Cause |"
    echo "|---:|---|"

    if [[ "${NB_MISSED}" -eq 0 ]]; then
        echo "| 0 | Aucune cause missed |"
    else
        awk -F'\t' '
        NR <= 10 {
            gsub(/\|/, "\\|", $2);
            print "| " $1 " | " $2 " |";
        }
        ' "${TOP_MISSED}"
    fi
} > "${MD_REPORT}"

echo "Rapport genere: ${MD_REPORT}"
echo "CSV genere: ${CSV_REPORT}"
echo "CSV prioritaire genere: ${CSV_FOCUS}"
echo
echo "Apercu:"
sed -n '1,40p' "${MD_REPORT}"

rm -rf "${TMP_DIR}"
