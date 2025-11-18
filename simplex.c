#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

// ========================
//      CONSTANTES
// ========================
#define EPS 1e-9
#define MAX_ITERS 10000

// ========================
//  WIDGETS
// ========================
GtkWidget *main_window;
GtkWidget *box_model;
GtkWidget *spin_vars;
GtkWidget *spin_constraints;
GtkWidget *second_window;
GtkWidget *check_intermediate;

GPtrArray *entry_var_names;
GPtrArray *labels_obj_vars;
GPtrArray *labels_con_vars;
GtkWidget *label_nonneg;
GtkWidget *combo_opt_global;

GPtrArray *entry_obj_coefs;
GPtrArray *entry_con_coefs;
GPtrArray *entry_rhs_list;
GPtrArray *combo_signs;

void on_generate_clicked(GtkButton *button, gpointer user_data);
void simplex(GtkWidget *btn, gpointer data);
void write_file(int vars, int cons, GtkWidget *box_model, int show_tables, const char *problem_name);
static double parse_coef(const char *text);

// ========================
// === Estructuras
// ========================
typedef struct
{
    int rows;         // 1 (fila objetivo) + cons
    int cols;         // 1 (Z) + vars + cons (holguras) + 1 (RHS)
    double **T;       // matriz [rows][cols]
    int *basis;       // índice de columna básica por cada fila (excepto fila 0). basis[r] = columna básica de la fila r (r>=1)
    int entering_col; // última col que entró
    int leaving_row;  // última fila que salió
    double BIG_M;            // valor numérico para computar (p.ej., 1e6)
    char **col_names;        // nombre por columna: "Z", "X_i", "s_i", "e_i", "a_i", "RHS"
    unsigned char *is_slack; // 1 si slack
    unsigned char *is_sur;   // 1 si exceso
    unsigned char *is_art;   // 1 si artificial

    // Mostrar M en tabla inicial (solo estética)
    int show_M_in_initial_row0;  // 1 si se quiere imprimir "M" en fila 0
    int *row0_M_sign;            // por columna: -1, 0, +1 (signo de M visible en fila 0)
    int row0_M_rhs_sign;         // por si quisieramos mostrar M en RHS (opc.)
} Tableau;

typedef struct
{
    Tableau *pre;        // snapshot ANTES de canonizar (para mostrar candidatos/razones)
    Tableau *post;       // snapshot DESPUÉS de canonizar (resultado de la iteración)
    int degenerate;      // 1 si en esta iteración se detectó degeneración (empate o razón 0)
    double min_ratio;    // razón mínima encontrada en la iteración (si aplica)
    GPtrArray *deg_rows; // array de rows (GINT_TO_POINTER(row)) que empataron o tienen ratio 0
} IterStep;

typedef struct
{
    GPtrArray *steps;        // Tableau* de cada iteración (se guardan copias)
    Tableau *initial;        // copia de la tabla inicial
    Tableau *final;          // copia de la tabla final (solución 1)
    Tableau *final_alt;      // copia de la tabla final alternativa (solución 2) si existe
    int status;              // 0=Óptimo, 1=No acotado, 2=No factible (no usamos en esta versión), 3=Degenerado (bandera informativa)
    int has_multiple;        // 1 si hay soluciones múltiples
    int degenerate;          // 1 si hubo razón mínima 0 o empates (informativo)
    double *solution;        // valores de variables de decisión (longitud = vars) - solución 1
    double *solution_alt;    // solución 2 (si existe)
    double **conv_solutions; // soluciones convexas adicionales (conv_count x vars)
    int conv_count;          // número de soluciones convexas (0 o 3)
    double z_value;          // valor óptimo de Z (para solución 1)
} SimplexRun;

// ========================
// === Utilidades Simplex
// ========================
static double **alloc_matrix(int r, int c)
{
    double **m = (double **)malloc(sizeof(double *) * r);
    for (int i = 0; i < r; i++)
    {
        m[i] = (double *)calloc(c, sizeof(double));
    }
    return m;
}
static void free_matrix(double **m, int r)
{
    if (!m)
        return;
    for (int i = 0; i < r; i++)
        free(m[i]);
    free(m);
}

static Tableau *tableau_new(int rows, int cols)
{
    Tableau *tb = (Tableau *)calloc(1, sizeof(Tableau));
    tb->rows = rows;
    tb->cols = cols;
    tb->BIG_M = 0.0;
    tb->col_names = NULL;
    tb->is_slack = tb->is_sur = tb->is_art = NULL;
    tb->show_M_in_initial_row0 = 0;
    tb->row0_M_sign = NULL;
    tb->row0_M_rhs_sign = 0;
    tb->T = alloc_matrix(rows, cols);
    tb->basis = (int *)malloc(sizeof(int) * rows);
    for (int i = 0; i < rows; i++)
        tb->basis[i] = -1;
    tb->entering_col = -1;
    tb->leaving_row = -1;
    return tb;
}


static Tableau *tableau_copy(const Tableau *src)
{
    if (!src) return NULL;

    // Crear nueva tabla con mismas dimensiones
    Tableau *t = tableau_new(src->rows, src->cols);

    // Copiar matriz numérica
    for (int i = 0; i < src->rows; i++) {
        memcpy(t->T[i], src->T[i], sizeof(double) * src->cols);
    }

    // Copiar base
    memcpy(t->basis, src->basis, sizeof(int) * src->rows);

    // Copiar campos escalares
    t->rows         = src->rows;
    t->cols         = src->cols;
    t->entering_col = src->entering_col;
    t->leaving_row  = src->leaving_row;
    t->BIG_M        = src->BIG_M;

    // === Copiar nombres de columnas ===
    if (src->col_names) {
        t->col_names = (char **)calloc(src->cols, sizeof(char *));
        for (int j = 0; j < src->cols; j++) {
            if (src->col_names[j]) {
                t->col_names[j] = strdup(src->col_names[j]);
            }
        }
    }

    // === Copiar flags de tipo de columna ===
    if (src->is_slack) {
        t->is_slack = (unsigned char *)calloc(src->cols, sizeof(unsigned char));
        memcpy(t->is_slack, src->is_slack, src->cols * sizeof(unsigned char));
    }

    if (src->is_sur) {
        t->is_sur = (unsigned char *)calloc(src->cols, sizeof(unsigned char));
        memcpy(t->is_sur, src->is_sur, src->cols * sizeof(unsigned char));
    }

    if (src->is_art) {
        t->is_art = (unsigned char *)calloc(src->cols, sizeof(unsigned char));
        memcpy(t->is_art, src->is_art, src->cols * sizeof(unsigned char));
    }

    // === Por ahora NO queremos mostrar M en las copias ===
    // (evitamos la lógica estética de M hasta que la queramos arreglar bien)
    t->show_M_in_initial_row0 = 0;
    t->row0_M_sign            = NULL;
    t->row0_M_rhs_sign        = 0;

    return t;
}

// Imprime un nombre de columna de forma segura en LaTeX
// - Si ya contiene '$', se asume que es código LaTeX y se imprime tal cual.
// - Si no contiene '$', se envuelve en $ ... $ y se escapan solo caracteres realmente peligrosos.
static void latex_escape_and_print(FILE *f, const char *s)
{
    // Si el nombre ya contiene $, asumimos que el usuario escribió LaTeX (ej: "$X_1$")
    if (strchr(s, '$') != NULL) {
        fputs(s, f);
        return;
    }

    // Si no tiene $, lo ponemos en modo matemático nosotros
    fputc('$', f);
    for (; *s; s++) {
        char c = *s;
        if (c == '%') {
            fprintf(f, "\\%%");
        } else if (c == '&') {
            fprintf(f, "\\&");
        } else if (c == '#') {
            fprintf(f, "\\#");
        } else if (c == '$') {
            fprintf(f, "\\$");
        } else {
            // OJO: aquí NO escapamos '_' porque estamos en modo matemático y queremos subíndices
            fputc(c, f);
        }
    }
    fputc('$', f);
}



static void tableau_free(Tableau *tb)
{
    if (!tb)
        return;
    if (tb->col_names){ for(int j=0;j<tb->cols;j++) free(tb->col_names[j]); free(tb->col_names); }
    free(tb->is_slack);
    free(tb->is_sur);
    free(tb->is_art);
    free(tb->row0_M_sign);
    free_matrix(tb->T, tb->rows);
    free(tb->basis);
    free(tb);
}

static IterStep *iterstep_new(Tableau *pre, Tableau *post)
{
    IterStep *s = (IterStep *)calloc(1, sizeof(IterStep));
    s->pre = pre;
    s->post = post;
    s->degenerate = 0;
    s->min_ratio = 0.0;
    s->deg_rows = NULL;
    return s;
}
static void iterstep_free(IterStep *s)
{
    if (!s)
        return;
    if (s->pre)
        tableau_free(s->pre);
    if (s->post)
        tableau_free(s->post);
    if (s->deg_rows)
        g_ptr_array_free(s->deg_rows, FALSE); // elementos son GINT_TO_POINTER, no free individual
    free(s);
}

static SimplexRun *simplex_run_new(void)
{
    SimplexRun *r = (SimplexRun *)calloc(1, sizeof(SimplexRun));
    r->steps = g_ptr_array_new_with_free_func((GDestroyNotify)iterstep_free);
    r->status = 0;
    r->has_multiple = 0;
    r->degenerate = 0;
    r->solution = NULL;
    r->solution_alt = NULL;
    r->conv_solutions = NULL;
    r->conv_count = 0;
    r->z_value = 0.0;
    r->final = NULL;
    r->final_alt = NULL;
    return r;
}
static void simplex_run_free(SimplexRun *run)
{
    if (!run)
        return;
    if (run->initial)
        tableau_free(run->initial);
    if (run->final)
        tableau_free(run->final);
    if (run->final_alt)
        tableau_free(run->final_alt);
    if (run->steps)
        g_ptr_array_free(run->steps, TRUE); // iterstep_free liberará tablas
    if (run->solution)
        free(run->solution);
    if (run->solution_alt)
        free(run->solution_alt);
    if (run->conv_solutions)
    {
        for (int i = 0; i < run->conv_count; i++)
            free(run->conv_solutions[i]);
        free(run->conv_solutions);
    }
    free(run);
}

// Helper para detectar degeneración y empates en razón mínima
static GPtrArray *collect_degeneracy_info(const Tableau *tb, int entering_col, double *out_min_ratio)
{
    GPtrArray *arr = g_ptr_array_new();
    double min_ratio = 1e300;
    // primera pasada: encontrar razón mínima entre a>EPS
    for (int i = 1; i < tb->rows; i++)
    {
        double a = tb->T[i][entering_col];
        if (a > EPS)
        {
            double rhs = tb->T[i][tb->cols - 1];
            double r = rhs / a;
            if (r < min_ratio)
                min_ratio = r;
        }
    }
    if (min_ratio == 1e300)
    {
        if (out_min_ratio)
            *out_min_ratio = 0.0;
        return arr; // vacío -> sin filas (posible no acotado)
    }
    // segunda pasada: recoger filas con ratio ~= min_ratio o ratio <= EPS (degeneración)
    for (int i = 1; i < tb->rows; i++)
    {
        double a = tb->T[i][entering_col];
        if (a > EPS)
        {
            double rhs = tb->T[i][tb->cols - 1];
            double r = rhs / a;
            if (fabs(r - min_ratio) <= EPS || r <= EPS)
            {
                g_ptr_array_add(arr, GINT_TO_POINTER(i));
            }
        }
    }
    if (out_min_ratio)
        *out_min_ratio = min_ratio;
    return arr;
}

// ========================
// === Construcción de la tabla inicial (Gran M)
// ========================
static Tableau *build_initial_tableau_from_inputs(int vars, int cons, int *unsupported_sign_found, int opt_type)
{
    // 1) Primera pasada: contar slacks, excesos, artificiales
    int cnt_slack = 0, cnt_sur = 0, cnt_art = 0;
    for (int i = 0; i < cons; i++) {
        GtkWidget *combo = g_ptr_array_index(combo_signs, i);
        int sign = gtk_combo_box_get_active(GTK_COMBO_BOX(combo)); // 0<=, 1=, 2>=
        if (sign == 0) cnt_slack++;
        else if (sign == 1) cnt_art++;          // '=' -> artificial
        else if (sign == 2) { cnt_sur++; cnt_art++; } // '>=' -> exceso + artificial
    }

    // 2) Definir dimensiones con todas las columnas
    int cols = 1 /*Z*/ + vars + cnt_slack + cnt_sur + cnt_art + 1 /*RHS*/;
    int rows = cons + 1;

    Tableau *tb = tableau_new(rows, cols);

    // Guarda M (numérico para cómputo) y meta
    tb->BIG_M = 1e6; // puedes exponerlo en la UI si deseas
    tb->col_names = (char**)calloc(cols, sizeof(char*));
    tb->is_slack = (unsigned char*)calloc(cols, 1);
    tb->is_sur   = (unsigned char*)calloc(cols, 1);
    tb->is_art   = (unsigned char*)calloc(cols, 1);
    tb->row0_M_sign = (int*)calloc(cols, sizeof(int));
    tb->show_M_in_initial_row0 = 1; // queremos mostrar M en la tabla inicial

    // 3) Nombres de columnas
    int cZ = 0; tb->col_names[cZ] = strdup("Z");
    int c = 1;
    for (int j = 0; j < vars; j++) {
        const char *vn = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, j)));
        tb->col_names[c++] = strdup(vn && *vn ? vn : "x");
    }
    int start_slack = c;
    for (int s = 0; s < cnt_slack; s++) {
        char buf[32]; snprintf(buf, sizeof(buf), "$s_{%d}$", s+1);

        tb->col_names[c] = strdup(buf); tb->is_slack[c] = 1; c++;
    }
    int start_sur = c;
    for (int e = 0; e < cnt_sur; e++) {
        char buf[16]; snprintf(buf, sizeof(buf), "$e_{%d}$", e+1);
        tb->col_names[c] = strdup(buf); tb->is_sur[c] = 1; c++;
    }
    int start_art = c;
    for (int a = 0; a < cnt_art; a++) {
        char buf[16]; snprintf(buf, sizeof(buf), "$a_{%d}$", a+1);
        tb->col_names[c] = strdup(buf); tb->is_art[c] = 1; c++;
    }
    int cRHS = cols - 1; tb->col_names[cRHS] = strdup("RHS");

    // 4) Fila 0: Z y FO
    tb->T[0][0] = 1.0;
    for (int j = 0; j < vars; j++) {
        const char *coef_str = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_obj_coefs, j)));
        double cj = parse_coef(coef_str);
        tb->T[0][1 + j] = -cj; // como ya usabas (Max usa "más negativo", Min "más positivo")
    }

    // 5) Llenar restricciones
    *unsupported_sign_found = 0;
    int idx_coef = 0;
    int pos_slack = 0, pos_sur = 0, pos_art = 0;

    for (int i = 0; i < cons; i++) {
        // Coeficientes de x
        for (int j = 0; j < vars; j++) {
            const char *a_str = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_con_coefs, idx_coef++)));
            tb->T[i+1][1 + j] = parse_coef(a_str);
        }
        // RHS
        GtkWidget *rhs = g_ptr_array_index(entry_rhs_list, i);
        tb->T[i+1][cRHS] = parse_coef(gtk_entry_get_text(GTK_ENTRY(rhs)));

        // Signo
        GtkWidget *combo = g_ptr_array_index(combo_signs, i);
        int sign = gtk_combo_box_get_active(GTK_COMBO_BOX(combo)); // 0<=,1=,2>=

        if (sign == 0) {
            // <= : slack
            int col_s = start_slack + (pos_slack++);
            tb->T[i+1][col_s] = 1.0;
            tb->basis[i+1] = col_s; // base inicial
        } else if (sign == 1) {
            // = : artificial
            int col_a = start_art + (pos_art++);
            tb->T[i+1][col_a] = 1.0;
            tb->basis[i+1] = col_a;
            // penalización M en fila 0
            tb->row0_M_sign[col_a] = (opt_type==0 ? +1 : -1); // Max:+M, Min:-M
            tb->T[0][col_a] += (opt_type==0 ? +tb->BIG_M : -tb->BIG_M);
        } else if (sign == 2) {
            // >= : exceso -1 y artificial +1
            int col_e = start_sur + (pos_sur++);
            int col_a = start_art + (pos_art++);
            tb->T[i+1][col_e] = -1.0;
            tb->T[i+1][col_a] = 1.0;
            tb->basis[i+1] = col_a;
            // penalización M en fila 0
            tb->row0_M_sign[col_a] = (opt_type==0 ? +1 : -1);
            tb->T[0][col_a] += (opt_type==0 ? +tb->BIG_M : -tb->BIG_M);
        } else {
            *unsupported_sign_found = 1;
        }
    }

    // 6) Canonizar fila 0 respecto a variables BÁSICAS artificiales
    //    (si una artificial está en la base y fila0 tiene coef ≠ 0 en esa col: fila0 -= coef * fila_básica)
    for (int r = 1; r < rows; r++) {
        int bc = tb->basis[r];
        if (bc >= 0 && tb->is_art[bc]) {
            double c0 = tb->T[0][bc];
            if (fabs(c0) > EPS) {
                for (int j = 0; j < cols; j++) tb->T[0][j] -= c0 * tb->T[r][j];
                // tras canonizar, el símbolo "M" ya no aparece algebraicamente en la numérica,
                // pero seguimos mostrando "M" en la tabla INICIAL a modo visual:
            }
        }
    }

    return tb;
}


// ========================
// Selección de columna que entra
//   max: más NEGATIVO en fila 0 (entre x y s no básicos, normalmente x)
//   min: más POSITIVO en fila 0
// Retorna -1 si ya óptimo según el criterio
// ========================
static int choose_entering_col(const Tableau *tb, int vars, int cons, int is_max)
{
    int start = 1;                             // después de Z
    int end = 1 + vars + cons - 1;             // última col de slack; RHS es la última col
    double best_val = is_max ? 1e100 : -1e100; // max: buscamos el más negativo => inicial alto; min: más positivo => inicial bajo
    int best_col = -1;

    for (int j = start; j < tb->cols - 1; j++)
    {
        double v = tb->T[0][j];
        if (is_max)
        {
            if (v < best_val - EPS)
            { // más negativo
                best_val = v;
                best_col = j;
            }
        }
        else
        { // min
            if (v > best_val + EPS)
            { // más positivo
                best_val = v;
                best_col = j;
            }
        }
    }

    // Criterio de optimalidad:
    // max: si NO hay negativos => óptimo
    // min: si NO hay positivos => óptimo
    if (is_max)
    {
        if (best_col == -1 || best_val >= -EPS)
            return -1;
    }
    else
    {
        if (best_col == -1 || best_val <= EPS)
            return -1;
    }
    return best_col;
}

// ========================
// NUEVO: Selección de fila que sale (test de razón mínima)
// Solo razones positivas. Marca degeneración si hay 0 o empates cercanos
// ========================
static int choose_leaving_row(const Tableau *tb, int entering_col, int *degenerate_flag)
{
    int best_row = -1;
    double best_ratio = 0.0;
    *degenerate_flag = 0;

    for (int i = 1; i < tb->rows; i++)
    {
        double a = tb->T[i][entering_col];
        if (a > EPS)
        { // solo positivos
            double rhs = tb->T[i][tb->cols - 1];
            double ratio = rhs / a;
            if (best_row == -1 || ratio < best_ratio - EPS)
            {
                best_ratio = ratio;
                best_row = i;
            }
            else if (fabs(ratio - best_ratio) <= EPS)
            {
                *degenerate_flag = 1; // empate de fracciones
            }
            if (ratio <= EPS)
                *degenerate_flag = 1; // degeneración (sale con 0)
        }
    }
    return best_row; // -1 => No acotado
}

// ========================
// Pivoteo (Gauss-Jordan en la tabla)
// ========================
static void pivot(Tableau *tb, int pr, int pc)
{
    double piv = tb->T[pr][pc];
    // Normalizar fila pivote
    for (int j = 0; j < tb->cols; j++)
        tb->T[pr][j] /= piv;

    // Anular columna pivote en las demás filas
    for (int i = 0; i < tb->rows; i++)
    {
        if (i == pr)
            continue;
        double factor = tb->T[i][pc];
        if (fabs(factor) > EPS)
        {
            for (int j = 0; j < tb->cols; j++)
            {
                tb->T[i][j] -= factor * tb->T[pr][j];
            }
        }
    }
    tb->basis[pr] = pc;
    tb->entering_col = pc;
    tb->leaving_row = pr;
}

// ========================
// Detección de soluciones múltiples
// max: en óptimo, si existe costo reducido 0 en columna NO básica (x), hay múltiples
// min: análogo (fila 0 con 0 en col no básica)
// ========================
static int detect_multiple_solutions(const Tableau *tb)
{
    // revisamos columnas no básicas en fila 0
    for (int j = 1; j < tb->cols - 1; j++)
    {
        // si la columna no está en la base:
        int in_basis = 0;
        for (int i = 1; i < tb->rows; i++)
            if (tb->basis[i] == j)
            {
                in_basis = 1;
                break;
            }
        if (!in_basis)
        {
            if (fabs(tb->T[0][j]) <= EPS)
                return 1;
        }
    }
    return 0;
}

// ========================
// Ejecutar Símplex (max/min)
// Guarda pasos si show_steps != 0
// ========================
static SimplexRun *simplex_solve(Tableau *init, int vars, int cons, int is_max, int show_steps)
{
    SimplexRun *run = simplex_run_new();
    run->initial = tableau_copy(init);
    Tableau *tb = tableau_copy(init);

    int iters = 0;
    while (iters++ < MAX_ITERS)
    {
        int entering = choose_entering_col(tb, vars, cons, is_max);
        if (entering == -1)
        { // óptimo
            run->status = 0;
            break;
        }
        int deg_flag = 0;
        int leaving = choose_leaving_row(tb, entering, &deg_flag);
        if (leaving == -1)
        { // no acotado
            run->status = 1;
            tb->entering_col = entering;
            tb->leaving_row = -1;
            if (show_steps)
            {
                Tableau *pre_snap = tableau_copy(tb);
                pre_snap->entering_col = entering;
                pre_snap->leaving_row = -1;
                IterStep *it = iterstep_new(pre_snap, NULL);
                g_ptr_array_add(run->steps, it);
            }
            break;
        }
        if (deg_flag)
            run->degenerate = 1;

        /* guardar snapshot ANTES del pivote y luego el resultado como UNA iteración */
        if (show_steps)
        {
            Tableau *pre = tableau_copy(tb);
            pre->entering_col = entering;
            pre->leaving_row = leaving;

            GPtrArray *degrows = NULL;
            double min_ratio = 0.0;
            if (deg_flag)
            {
                degrows = collect_degeneracy_info(pre, entering, &min_ratio);
            }

            // aplicamos pivote sobre la tabla real
            pivot(tb, leaving, entering);

            Tableau *post = tableau_copy(tb);
            // post ya tiene entering/leaving asignados por pivot

            IterStep *it = iterstep_new(pre, post);
            if (deg_flag)
            {
                it->degenerate = 1;
                it->min_ratio = min_ratio;
                it->deg_rows = degrows;
            }
            g_ptr_array_add(run->steps, it);
        }
        else
        {
            // si no mostramos pasos, igual hacemos el pivote normal
            pivot(tb, leaving, entering);
        }
    }

    run->final = tableau_copy(tb);

    /* solución vectorial para final (solución 1) */
    run->solution = (double *)calloc(vars, sizeof(double));
    for (int i = 1; i < run->final->rows; i++)
    {
        int bc = run->final->basis[i];
        if (bc >= 1 && bc <= vars)
        {
            run->solution[bc - 1] = run->final->T[i][run->final->cols - 1];
        }
    }
    run->z_value = run->final->T[0][run->final->cols - 1];

    /* Detectar múltiples soluciones sobre la tabla final */
    run->has_multiple = detect_multiple_solutions(run->final);
    if (run->has_multiple)
    {
        /* Intentar pivotear una vez más en alguna columna no básica con costo reducido ~ 0 */
        int pivoted = 0;
        for (int j = 1; j < run->final->cols - 1 && !pivoted; j++)
        {
            int in_basis = 0;
            for (int ii = 1; ii < run->final->rows; ii++)
                if (run->final->basis[ii] == j)
                {
                    in_basis = 1;
                    break;
                }
            if (!in_basis && fabs(run->final->T[0][j]) <= EPS)
            {
                int dummy_deg = 0;
                int lv = choose_leaving_row(run->final, j, &dummy_deg);
                if (lv != -1)
                {
                    /* no modificamos run->final: creamos copia y pivotamos para obtener solución 2 */
                    Tableau *alt = tableau_copy(run->final);
                    pivot(alt, lv, j);
                    run->final_alt = alt;

                    /* solución alternativa */
                    run->solution_alt = (double *)calloc(vars, sizeof(double));
                    for (int ii = 1; ii < alt->rows; ii++)
                    {
                        int bc2 = alt->basis[ii];
                        if (bc2 >= 1 && bc2 <= vars)
                        {
                            run->solution_alt[bc2 - 1] = alt->T[ii][alt->cols - 1];
                        }
                    }

                    /* Construir 3 soluciones convexas: alfa * sol1 + (1 - alfa) * sol2 */
                    double alphas[3] = {0.25, 0.5, 0.75};
                    run->conv_count = 3;
                    run->conv_solutions = (double **)malloc(sizeof(double *) * run->conv_count);
                    for (int a = 0; a < run->conv_count; a++)
                    {
                        run->conv_solutions[a] = (double *)calloc(vars, sizeof(double));
                        double alpha = alphas[a];
                        for (int v = 0; v < vars; v++)
                        {
                            double s1 = run->solution ? run->solution[v] : 0.0;
                            double s2 = run->solution_alt ? run->solution_alt[v] : 0.0;
                            run->conv_solutions[a][v] = alpha * s1 + (1.0 - alpha) * s2;
                        }
                    }

                    pivoted = 1;
                }
            }
        }
        /* marca degeneración informativa si no se pudo pivotear adicionalmente */
        if (!pivoted)
            run->degenerate = 1;
    }

    tableau_free(tb);
    return run;
}

// ========================
// Impresores LaTeX
// ========================
static void latex_print_table(FILE *f, const Tableau *tb, int vars, int cons,
                              const char **var_names, int highlight_col, int highlight_row,
                              int show_ratio, int ratio_col)
{
    // --- Determinar visibilidad de columnas ---
    // Visible si:
    //  - Z y RHS siempre.
    //  - Variables de decisión y slacks y excesos siempre.
    //  - Artificiales SOLO si están en la base actual.
    int *visible = (int*)calloc(tb->cols, sizeof(int));
    for (int j = 0; j < tb->cols; j++) {
        if (j == 0 || j == tb->cols-1) { visible[j] = 1; continue; }
        if (tb->is_art && tb->is_art[j]) {
            int in_basis = 0;
            for (int i = 1; i < tb->rows; i++) if (tb->basis[i] == j) { in_basis = 1; break; }
            visible[j] = in_basis;
        } else {
            visible[j] = 1;
        }
    }

    // --- Razones (si aplica) ---
    if (show_ratio && ratio_col >= 0) {
        fprintf(f, "\\\\[0.2cm]\n");
        fprintf(f, "Cálculo de razones de valores positivos de la columna %d: \\\\\n\n", ratio_col);

        fprintf(f, "\\begin{tabular}{lr}\\toprule Fila & Razón \\\\\\midrule\n");
        for (int i = 1; i < tb->rows; i++) {
            double a = tb->T[i][ratio_col];
            if (a > EPS) {
                double rhs = tb->T[i][tb->cols - 1];
                double r = rhs / a;
                fprintf(f, "%d & %.6g \\\\\n", i, r);
            }
        }
        fprintf(f, "\\bottomrule\\end{tabular}\n\n\\vspace{0.4cm}\n");

        if (highlight_row >= 0) {
            double piv_val = tb->T[highlight_row][ratio_col];
            fprintf(f, "\\noindent Con %.6g como pivote en columna %d.\\\\\n", piv_val, ratio_col);
        } else {
            fprintf(f, "\\noindent No se encuentran candidatos a pivote en columna %d.\\\\\n", ratio_col);
        }
    }

    // --- Tabla ---
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\rowcolors{2}{white}{gray!10}\n");

    // contar visibles
    int vcount = 0;
    for (int j = 0; j < tb->cols; j++) if (visible[j]) vcount++;

    fprintf(f, "\\begin{tabular}{");
    for (int j = 0; j < vcount; j++) fprintf(f, "r");
    fprintf(f, "}\n\\toprule\n");

    // Encabezados
    int printed = 0;
    for (int j = 0; j < tb->cols; j++) {
        if (!visible[j]) continue;
        if (printed) fprintf(f, " & ");

        const char *name = (tb->col_names && tb->col_names[j])
                        ? tb->col_names[j]
                        : (j==0 ? "Z" : (j==tb->cols-1 ? "RHS" : "col"));

        if (j == highlight_col) {
            fprintf(f, "\\cellcolor{yellow!40}\\textbf{");
            latex_escape_and_print(f, name);
            fprintf(f, "}");
        } else {
            fprintf(f, "\\textbf{");
            latex_escape_and_print(f, name);
            fprintf(f, "}");
        }
        printed++;
    }
    fprintf(f, " \\\\\\midrule\n");

    // Filas
    for (int i = 0; i < tb->rows; i++) {
        printed = 0;
        for (int j = 0; j < tb->cols; j++) {
            if (!visible[j]) continue;

            if (printed) fprintf(f, " & ");
            int is_pivot = (i == highlight_row && j == highlight_col);

            // Mostrar "M" en fila 0 inicial (solo estética)
            if (i == 0 && tb->show_M_in_initial_row0 && tb->row0_M_sign && tb->row0_M_sign[j] != 0) {
                // Imprime algo tipo:  M   o   -M   o   (num ± M)
                double num = tb->T[i][j] - tb->row0_M_sign[j]*tb->BIG_M; // separar parte M
                if (is_pivot) fprintf(f, "\\cellcolor{yellow!40}");
                if (fabs(num) > 1e-9)
                    fprintf(f, "%.6g %c M", num, tb->row0_M_sign[j] > 0 ? '+' : '-');
                else
                    fprintf(f, "%sM", tb->row0_M_sign[j] > 0 ? "" : "-");
            } else {
                if (is_pivot) fprintf(f, "\\cellcolor{yellow!40}");
                fprintf(f, "%.6g", tb->T[i][j]);
            }

            printed++;
        }
        fprintf(f, "\\\\\n");
    }

    fprintf(f, "\\bottomrule\n\\end{tabular}\n");
    fprintf(f, "\\end{center}\n");

    free(visible);
}


static void latex_print_all_variables(FILE *f, const Tableau *tb,
                                      const char **var_names, int vars, int cons)
{
    fprintf(f, "\\section*{Valores de todas las variables}\n");
    fprintf(f, "\\begin{tabular}{lr}\\toprule Variable & Valor \\\\ \\midrule\n");

    // Recorremos columnas 1..vars (variables de decisión)
    for (int j = 1; j <= vars; j++)
    {
        double val = 0.0;
        for (int i = 1; i < tb->rows; i++)
        {
            if (tb->basis[i] == j)
            {
                val = tb->T[i][tb->cols - 1];
                break;
            }
        }
        fprintf(f, "%s & %.6g \\\\\n", var_names[j - 1], val);
    }

    // Recorremos holguras s_i (columnas vars+1 ... vars+cons)
    for (int j = vars + 1; j <= vars + cons; j++)
    {
        double val = 0.0;
        for (int i = 1; i < tb->rows; i++)
        {
            if (tb->basis[i] == j)
            {
                val = tb->T[i][tb->cols - 1];
                break;
            }
        }
        fprintf(f, "s_{%d} & %.6g \\\\\n", j - vars, val);
    }

    fprintf(f, "\\bottomrule\\end{tabular}\\\\[0.3cm]\n");
}

// ========================
//  HANDLERS GUI
// ========================

// --- Cerrar todo el programa ---
void on_close_main_window(GtkWidget *widget, gpointer user_data)
{
    gtk_main_quit();
}

// --- Cuando se cierra la segunda ventana ---
void on_second_window_destroy(GtkWidget *widget, gpointer user_data)
{
    gtk_widget_hide(second_window);
    gtk_widget_show_all(main_window);
}

// --- Botón "Cancelar" ---
void on_cancel_clicked(GtkButton *button, gpointer user_data)
{
    gtk_widget_hide(second_window);
    gtk_widget_show_all(main_window);
}

// --- Actualizar "Y con X1, X2 ≥ 0" ---
void update_nonneg_label()
{
    GString *vars_text = g_string_new("Y con ");
    for (int i = 0; i < entry_var_names->len; i++)
    {
        const char *name = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, i)));
        g_string_append(vars_text, name);
        if (i < entry_var_names->len - 1)
            g_string_append(vars_text, ", ");
    }
    g_string_append(vars_text, " ≥ 0");
    gtk_label_set_text(GTK_LABEL(label_nonneg), vars_text->str);
    g_string_free(vars_text, TRUE);
}

// --- Actualizar nombres en la segunda ventana ---
void on_var_name_changed(GtkEditable *editable, gpointer user_data)
{
    int index = GPOINTER_TO_INT(user_data);
    const char *new_name = gtk_entry_get_text(GTK_ENTRY(editable));

    if (index < labels_obj_vars->len)
    {
        GtkWidget *label = g_ptr_array_index(labels_obj_vars, index);
        gtk_label_set_text(GTK_LABEL(label), new_name);
    }

    for (int i = 0; i < labels_con_vars->len; i++)
    {
        GPtrArray *fila = g_ptr_array_index(labels_con_vars, i);
        if (index < fila->len)
        {
            GtkWidget *label = g_ptr_array_index(fila, index);
            gtk_label_set_text(GTK_LABEL(label), new_name);
        }
    }

    update_nonneg_label();
}

// --- Guardar ---
void on_save_clicked(GtkButton *button, gpointer user_data)
{
    GtkWidget *dialog;
    gchar *filename = NULL;

    dialog = gtk_file_chooser_dialog_new(
        "Guardar modelo",
        NULL,
        GTK_FILE_CHOOSER_ACTION_SAVE,
        "_Guardar", GTK_RESPONSE_ACCEPT,
        "_Cancelar", GTK_RESPONSE_CANCEL,
        NULL);

    gtk_file_chooser_set_do_overwrite_confirmation(GTK_FILE_CHOOSER(dialog), TRUE);
    GtkFileFilter *filter = gtk_file_filter_new();
    gtk_file_filter_add_pattern(filter, "*.txt");
    gtk_file_filter_set_name(filter, "Archivos de texto (*.txt)");
    gtk_file_chooser_add_filter(GTK_FILE_CHOOSER(dialog), filter);

    const char *problem_name = g_object_get_data(G_OBJECT(button), "problem_name");
    if (problem_name && strlen(problem_name) > 0)
    {
        char file_name[256];
        snprintf(file_name, sizeof(file_name), "%s.txt", problem_name);
        gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(dialog), file_name);
    }

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_ACCEPT)
    {
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));
        if (!g_str_has_suffix(filename, ".txt"))
        {
            char *with_ext = g_strconcat(filename, ".txt", NULL);
            g_free(filename);
            filename = with_ext;
        }
    }
    else
    {
        gtk_widget_destroy(dialog);
        return;
    }
    gtk_widget_destroy(dialog);

    FILE *f = fopen(filename, "w");
    if (!f)
    {
        g_free(filename);
        return;
    }

    // Nombre del problema
    fprintf(f, "%s\n", problem_name ? problem_name : "SinNombre");

    // Guardar tipo de optimización (0 = Max, 1 = Min)
    int opt_type = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_opt_global));
    fprintf(f, "%d\n", opt_type);

    // Cantidad de variables y restricciones
    fprintf(f, "%u\n", entry_var_names->len);
    fprintf(f, "%u\n", labels_con_vars->len);

    // Nombres de variables
    for (int i = 0; i < entry_var_names->len; i++)
        fprintf(f, "%s\n", gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, i))));

    // Coeficientes de Z
    for (int i = 0; i < entry_obj_coefs->len; i++)
        fprintf(f, "%s\n", gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_obj_coefs, i))));

    // Restricciones
    int vars = entry_var_names->len;
    int cons = labels_con_vars->len;
    int idx = 0;

    for (int r = 0; r < cons; r++)
    {
        for (int v = 0; v < vars; v++)
        {
            GtkWidget *entry = g_ptr_array_index(entry_con_coefs, idx++);
            fprintf(f, "%s\n", gtk_entry_get_text(GTK_ENTRY(entry)));
        }
        GtkWidget *combo = g_ptr_array_index(combo_signs, r);
        GtkWidget *rhs = g_ptr_array_index(entry_rhs_list, r);
        fprintf(f, "%d\n", gtk_combo_box_get_active(GTK_COMBO_BOX(combo)));
        fprintf(f, "%s\n", gtk_entry_get_text(GTK_ENTRY(rhs)));
    }

    fclose(f);
    g_free(filename);
}

// --- Cargar ---
void on_load_clicked(GtkButton *button, gpointer user_data)
{
    GtkWidget *dialog;
    gchar *filename = NULL;

    dialog = gtk_file_chooser_dialog_new(
        "Cargar modelo",
        NULL,
        GTK_FILE_CHOOSER_ACTION_OPEN,
        "_Abrir", GTK_RESPONSE_OK,
        "_Cancelar", GTK_RESPONSE_CANCEL,
        NULL);

    if (gtk_dialog_run(GTK_DIALOG(dialog)) == GTK_RESPONSE_OK)
        filename = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(dialog));

    gtk_widget_destroy(dialog);
    if (!filename)
        return;

    FILE *f = fopen(filename, "r");
    if (!f)
    {
        g_free(filename);
        return;
    }

    char problem_name[128];
    fgets(problem_name, sizeof(problem_name), f);
    problem_name[strcspn(problem_name, "\r\n")] = 0;

    int opt_type;
    fscanf(f, "%d\n", &opt_type);

    int vars, cons;
    fscanf(f, "%d\n%d\n", &vars, &cons);

    GtkWidget *entry_problem_name = GTK_WIDGET(gtk_builder_get_object(user_data, "entry_problem_name"));
    gtk_entry_set_text(GTK_ENTRY(entry_problem_name), problem_name);

    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin_vars), vars);
    gtk_spin_button_set_value(GTK_SPIN_BUTTON(spin_constraints), cons);

    on_generate_clicked(NULL, user_data);

    gtk_combo_box_set_active(GTK_COMBO_BOX(combo_opt_global), opt_type);

    char buf[128];
    for (int i = 0; i < vars; i++)
    {
        fgets(buf, sizeof(buf), f);
        buf[strcspn(buf, "\r\n")] = 0;
        gtk_entry_set_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, i)), buf);
    }

    for (int i = 0; i < vars; i++)
    {
        fgets(buf, sizeof(buf), f);
        buf[strcspn(buf, "\r\n")] = 0;
        gtk_entry_set_text(GTK_ENTRY(g_ptr_array_index(entry_obj_coefs, i)), buf);
    }

    int idx = 0;
    for (int r = 0; r < cons; r++)
    {
        for (int v = 0; v < vars; v++)
        {
            fgets(buf, sizeof(buf), f);
            buf[strcspn(buf, "\r\n")] = 0;
            gtk_entry_set_text(GTK_ENTRY(g_ptr_array_index(entry_con_coefs, idx++)), buf);
        }
        int sign;
        fscanf(f, "%d\n", &sign);
        GtkWidget *combo = g_ptr_array_index(combo_signs, r);
        gtk_combo_box_set_active(GTK_COMBO_BOX(combo), sign);

        fgets(buf, sizeof(buf), f);
        buf[strcspn(buf, "\r\n")] = 0;
        GtkWidget *rhs = g_ptr_array_index(entry_rhs_list, r);
        gtk_entry_set_text(GTK_ENTRY(rhs), buf);
    }

    fclose(f);
    g_free(filename);
}

// Borrar cuando se implemente la gran M
void on_exec_clicked(GtkButton *button, gpointer user_data)
{

    // Si todo está bien → continuar con la ejecución original:
    simplex(GTK_WIDGET(button), user_data);
}

// --- Generar modelo ---
void on_generate_clicked(GtkButton *button, gpointer user_data)
{
    GtkWidget *entry_problem_name_widget = GTK_WIDGET(gtk_builder_get_object(user_data, "entry_problem_name"));
    const char *problem_name = gtk_entry_get_text(GTK_ENTRY(entry_problem_name_widget));

    if (!problem_name || strlen(problem_name) == 0)
    {
        GtkWidget *dialog = gtk_message_dialog_new(GTK_WINDOW(main_window),
                                                   GTK_DIALOG_MODAL,
                                                   GTK_MESSAGE_WARNING,
                                                   GTK_BUTTONS_OK,
                                                   "Por favor, ingrese un nombre para el problema.");
        gtk_dialog_run(GTK_DIALOG(dialog));
        gtk_widget_destroy(dialog);
        return;
    }

    int vars = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(spin_vars));
    int cons = gtk_spin_button_get_value_as_int(GTK_SPIN_BUTTON(spin_constraints));

    gtk_widget_hide(main_window);

    second_window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
    gtk_window_set_title(GTK_WINDOW(second_window), problem_name);
    gtk_window_set_default_size(GTK_WINDOW(second_window), 600, 400);
    g_signal_connect(second_window, "destroy", G_CALLBACK(on_second_window_destroy), NULL);

    GtkWidget *main_vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL, 10);
    gtk_container_add(GTK_CONTAINER(second_window), main_vbox);
    gtk_container_set_border_width(GTK_CONTAINER(main_vbox), 10);

    entry_var_names = g_ptr_array_new();
    entry_obj_coefs = g_ptr_array_new();
    entry_con_coefs = g_ptr_array_new();
    entry_rhs_list = g_ptr_array_new();
    combo_signs = g_ptr_array_new();
    labels_obj_vars = g_ptr_array_new();
    labels_con_vars = g_ptr_array_new();

    // Nombres de variables
    GtkWidget *label_names = gtk_label_new("Nombres de las variables:");
    gtk_widget_set_halign(label_names, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(main_vbox), label_names, FALSE, FALSE, 0);

    GtkWidget *hbox_names = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
    gtk_box_pack_start(GTK_BOX(main_vbox), hbox_names, FALSE, FALSE, 0);

    for (int i = 1; i <= vars; i++)
    {
        GtkWidget *entry = gtk_entry_new();
        char def_name[16];
        snprintf(def_name, sizeof(def_name), "$X_%d$", i);
        gtk_entry_set_text(GTK_ENTRY(entry), def_name);
        gtk_entry_set_width_chars(GTK_ENTRY(entry), 5);
        gtk_box_pack_start(GTK_BOX(hbox_names), entry, FALSE, FALSE, 5);
        g_ptr_array_add(entry_var_names, entry);
    }

    // Modelo
    GtkWidget *label_model = gtk_label_new("Modelo:");
    gtk_widget_set_halign(label_model, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(main_vbox), label_model, FALSE, FALSE, 0);

    GtkWidget *hbox_obj = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
    combo_opt_global = gtk_combo_box_text_new();
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo_opt_global), "Maximizar");
    gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo_opt_global), "Minimizar");
    gtk_combo_box_set_active(GTK_COMBO_BOX(combo_opt_global), 0);
    gtk_box_pack_start(GTK_BOX(hbox_obj), combo_opt_global, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(hbox_obj), gtk_label_new("Z ="), FALSE, FALSE, 5);

    for (int i = 0; i < vars; i++)
    {
        if (i > 0)
            gtk_box_pack_start(GTK_BOX(hbox_obj), gtk_label_new("+"), FALSE, FALSE, 3);

        GtkWidget *entry_coef = gtk_entry_new();
        gtk_entry_set_width_chars(GTK_ENTRY(entry_coef), 5);
        gtk_entry_set_text(GTK_ENTRY(entry_coef), "0");
        gtk_box_pack_start(GTK_BOX(hbox_obj), entry_coef, FALSE, FALSE, 5);
        g_ptr_array_add(entry_obj_coefs, entry_coef);

        const char *var_name = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, i)));
        GtkWidget *label_var = gtk_label_new(var_name);
        gtk_box_pack_start(GTK_BOX(hbox_obj), label_var, FALSE, FALSE, 5);
        g_ptr_array_add(labels_obj_vars, label_var);
    }

    gtk_box_pack_start(GTK_BOX(main_vbox), hbox_obj, FALSE, FALSE, 0);

    GtkWidget *label_sujeto = gtk_label_new("Sujeto a:");
    gtk_widget_set_halign(label_sujeto, GTK_ALIGN_START);
    gtk_box_pack_start(GTK_BOX(main_vbox), label_sujeto, FALSE, FALSE, 5);

    for (int r = 0; r < cons; r++)
    {
        GtkWidget *hbox_con = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
        GPtrArray *fila_labels = g_ptr_array_new();

        for (int v = 0; v < vars; v++)
        {
            if (v > 0)
                gtk_box_pack_start(GTK_BOX(hbox_con), gtk_label_new("+"), FALSE, FALSE, 3);

            GtkWidget *entry_coef = gtk_entry_new();
            gtk_entry_set_width_chars(GTK_ENTRY(entry_coef), 5);
            gtk_entry_set_text(GTK_ENTRY(entry_coef), "0");
            gtk_box_pack_start(GTK_BOX(hbox_con), entry_coef, FALSE, FALSE, 5);
            g_ptr_array_add(entry_con_coefs, entry_coef);

            const char *var_name = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, v)));
            GtkWidget *label_var = gtk_label_new(var_name);
            gtk_box_pack_start(GTK_BOX(hbox_con), label_var, FALSE, FALSE, 5);
            g_ptr_array_add(fila_labels, label_var);
        }

        GtkWidget *combo_sign = gtk_combo_box_text_new();
        gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo_sign), "<=");
        gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo_sign), "=");
        gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combo_sign), ">=");
        gtk_combo_box_set_active(GTK_COMBO_BOX(combo_sign), 0);
        gtk_box_pack_start(GTK_BOX(hbox_con), combo_sign, FALSE, FALSE, 5);
        g_ptr_array_add(combo_signs, combo_sign);

        GtkWidget *entry_rhs = gtk_entry_new();
        gtk_entry_set_width_chars(GTK_ENTRY(entry_rhs), 5);
        gtk_entry_set_text(GTK_ENTRY(entry_rhs), "0");
        gtk_box_pack_start(GTK_BOX(hbox_con), entry_rhs, FALSE, FALSE, 5);
        g_ptr_array_add(entry_rhs_list, entry_rhs);

        gtk_box_pack_start(GTK_BOX(main_vbox), hbox_con, FALSE, FALSE, 0);
        g_ptr_array_add(labels_con_vars, fila_labels);
    }

    label_nonneg = gtk_label_new("");
    gtk_box_pack_start(GTK_BOX(main_vbox), label_nonneg, FALSE, FALSE, 10);
    update_nonneg_label();

    // Checkbox "Mostrar tablas intermedias"
    check_intermediate = gtk_check_button_new_with_label("Mostrar tablas intermedias");
    gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check_intermediate), FALSE);
    gtk_box_pack_start(GTK_BOX(main_vbox), check_intermediate, FALSE, FALSE, 10);

    GtkWidget *hbox_buttons = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 10);
    GtkWidget *btn_cancel = gtk_button_new_with_label("Cancelar");
    GtkWidget *btn_save = gtk_button_new_with_label("Guardar");
    GtkWidget *btn_exec = gtk_button_new_with_label("Ejecutar");

    gtk_box_pack_start(GTK_BOX(hbox_buttons), btn_cancel, FALSE, FALSE, 5);
    gtk_box_pack_end(GTK_BOX(hbox_buttons), btn_exec, FALSE, FALSE, 5);
    gtk_box_pack_end(GTK_BOX(hbox_buttons), btn_save, FALSE, FALSE, 5);
    gtk_box_pack_start(GTK_BOX(main_vbox), hbox_buttons, FALSE, FALSE, 10);

    g_object_set_data_full(G_OBJECT(btn_save), "problem_name", g_strdup(problem_name), g_free);
    g_object_set_data_full(G_OBJECT(btn_exec), "problem_name", g_strdup(problem_name), g_free);
    g_signal_connect(btn_cancel, "clicked", G_CALLBACK(on_cancel_clicked), NULL);
    g_signal_connect(btn_save, "clicked", G_CALLBACK(on_save_clicked), NULL);
    GtkWidget **data = g_new(GtkWidget *, 3);
    data[0] = box_model;
    data[1] = GINT_TO_POINTER(vars);
    data[2] = GINT_TO_POINTER(cons);

    g_signal_connect_data(btn_exec, "clicked", G_CALLBACK(on_exec_clicked), data,
                          (GClosureNotify)g_free, 0);

    for (int i = 0; i < entry_var_names->len; i++)
    {
        GtkWidget *entry = g_ptr_array_index(entry_var_names, i);
        g_signal_connect(entry, "changed", G_CALLBACK(on_var_name_changed), GINT_TO_POINTER(i));
    }

    gtk_widget_show_all(second_window);
}

// obtener nombres de variables en arreglo C
static void collect_var_names(int vars, const char **names_buf)
{
    for (int i = 0; i < vars; i++)
    {
        names_buf[i] = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, i)));
    }
}

static double parse_coef(const char *text)
{
    if (!text)
        return 0.0;

    // Elimina espacios iniciales
    while (isspace((unsigned char)*text))
        text++;

    // Ver si es fracción
    const char *slash = strchr(text, '/');
    if (slash)
    {
        char num_buf[64], den_buf[64];

        size_t len_num = slash - text;
        size_t len_den = strlen(slash + 1);

        if (len_num >= sizeof(num_buf))
            len_num = sizeof(num_buf) - 1;
        if (len_den >= sizeof(den_buf))
            len_den = sizeof(den_buf) - 1;

        memcpy(num_buf, text, len_num);
        num_buf[len_num] = '\0';

        memcpy(den_buf, slash + 1, len_den);
        den_buf[len_den] = '\0';

        double num = atof(num_buf);
        double den = atof(den_buf);

        if (fabs(den) < 1e-12)
            return 0.0;

        return num / den;
    }

    // Si no hay "/", usar atof normal
    return atof(text);
}

static void print_term(FILE *f, const char *coef_str, const char *var, int is_first)
{
    double coef = parse_coef(coef_str);

    if (fabs(coef) < 1e-12)
        return; // coeficiente = 0 → NO imprimir

    // Signo
    if (!is_first)
    {
        if (coef > 0)
            fprintf(f, " + ");
        else
            fprintf(f, " - ");
    }
    else
    {
        if (coef < 0)
            fprintf(f, "-");
    }

    double abs_coef = fabs(coef);

    // Si coef = 1 → imprimir solo variable
    if (fabs(abs_coef - 1.0) < 1e-12)
        fprintf(f, "%s", var);
    else
        fprintf(f, "%.6g %s", abs_coef, var);
}

// --- Escritura del Latex ---
void write_file(int vars, int cons, GtkWidget *box_model, int show_tables, const char *problem_name)
{
    // 1) Construir tabla inicial del Símplex
    int unsupported_signs = 0;
    int opt_type = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_opt_global)); // 0=max, 1=min
    Tableau *initial = build_initial_tableau_from_inputs(vars, cons, &unsupported_signs, opt_type);

    // 2) Resolver
    SimplexRun *run = simplex_solve(initial, vars, cons, opt_type == 0, show_tables != 0);

    // 3) Abrir .tex
    FILE *f = fopen("simplex_output.tex", "w");
    if (!f)
    {
        perror("fopen tex");
        tableau_free(initial);
        simplex_run_free(run);
        return;
    }

    fprintf(f, "%% Archivo generado automáticamente por simplex.c\n");
    fprintf(f, "\\documentclass[12pt]{article}\n");
    fprintf(f, "\\usepackage[utf8]{inputenc}\n");
    fprintf(f, "\\usepackage{graphicx}\n");
    fprintf(f, "\\usepackage{array,booktabs}\n");
    fprintf(f, "\\usepackage[table]{xcolor}\n");
    fprintf(f, "\\usepackage{longtable}\n");
    fprintf(f, "\\usepackage{geometry}\n");
    fprintf(f, "\\usepackage{pdflscape}\n");
    fprintf(f, "\\usepackage{tikz}\n");
    fprintf(f, "\\usepackage{pgfplots}\n");
    fprintf(f, "\\pgfplotsset{compat=1.18}\n");
    fprintf(f, "\\geometry{margin=0.8in}\n");
    fprintf(f, "\\title{Proyecto 4: SIMPLEX}\n");
    fprintf(f, "\\author{Investigación de Operaciones}\n");
    fprintf(f, "\\date{}\n");
    fprintf(f, "\\begin{document}\n");

    // --- Portada ---
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "{\\large Instituto Tecnológico de Costa Rica\\\\[1cm]\n");
    fprintf(f, "\\includegraphics[width=0.4\\textwidth]{TEC.png}\\\\[2cm]\n");
    fprintf(f, "{\\LARGE \\textbf{Proyecto 4: SIMPLEX}}\\\\[2cm]\n");
    fprintf(f, "{\\large Investigación de Operaciones\\\\[2cm]\n");
    fprintf(f, "{\\large Profesor: }\\\\[1cm]\n");
    fprintf(f, "{\\large Francisco Jose Torres Rojas}\\\\[2cm]\n");
    fprintf(f, "{\\large Integrantes: }\\\\[1cm]\n");
    fprintf(f, "{\\large Jose Pablo Fernandez Jimenez - 2023117752}\\\\[1cm]\n");
    fprintf(f, "{\\large Diego Durán Rodríguez - 2022437509}\\\\[2cm]\n");
    fprintf(f, "{\\large Segundo semestre 2025\\\\[1cm]\n");
    fprintf(f, "\\end{center}\n\\newpage\n");

    // --- Algoritmo Simplex---
    fprintf(f, "\\section*{Algoritmo SIMPLEX}\n");
    fprintf(f,
            "El \\textbf{algoritmo SIMPLEX} es uno de los métodos más importantes y utilizados en el campo de la optimización lineal. "
            "Fue desarrollado por George Bernard Dantzig en 1947, en el contexto de investigaciones relacionadas con la planificación "
            "logística y de recursos durante la posguerra. Dantzig, un matemático y científico estadounidense, ideó este método como una "
            "herramienta para resolver problemas de programación lineal, un área que busca optimizar (maximizar o minimizar) una función "
            "objetivo sujeta a un conjunto de restricciones lineales. \n\n"

            "El surgimiento del algoritmo Simplex marcó un antes y un después en la optimización matemática. Antes de su creación, no existía "
            "un procedimiento general y sistemático que permitiera resolver eficientemente problemas de gran escala con múltiples variables y "
            "restricciones. Dantzig propuso un enfoque geométrico basado en la observación de que la solución óptima de un problema lineal se "
            "encuentra en uno de los vértices o puntos extremos del poliedro factible, es decir, del conjunto de soluciones que cumplen todas "
            "las restricciones del problema. \n\n"

            "El método Simplex avanza de un vértice a otro a través de las aristas del poliedro, mejorando progresivamente el valor de la "
            "función objetivo hasta encontrar el óptimo. Cada movimiento corresponde a un cambio de una variable básica en la solución, lo "
            "que permite al algoritmo recorrer el espacio factible de manera ordenada y eficiente. \n\n"

            "\\textbf{Entre las principales propiedades del algoritmo Simplex destacan las siguientes:}\n"
            "\\begin{itemize}\n"
            "  \\item \\textit{Eficiencia práctica:} aunque en teoría su complejidad puede ser exponencial en el peor de los casos, en la práctica el algoritmo es extremadamente eficiente y puede resolver problemas con miles de variables y restricciones en tiempos incluso lineales.\n"
            "  \\item \\textit{Interpretación geométrica clara:} El procedimiento del Simplex se basa en conceptos geométricos simples, lo que facilita su comprensión y visualización en espacios de baja dimensión.\n"
            "  \\item \\textit{Importancia histórica y teórica:} el Simplex no solo revolucionó la programación lineal, sino que también sentó las bases para la aparición de otros métodos de optimización, como los algoritmos de punto interior y técnicas modernas de optimización convexa.\n"
            "\\end{itemize}\n\n"

            "\\bigskip\n");

    // --- Problema original ---
    fprintf(f, "\\section*{Problema original}\n");
    fprintf(f, "Nombre del problema: \\textbf{%s}\\\\[0.3cm]\n", problem_name);
    fprintf(f, "El problema original se puede formular como un problema de programación lineal, donde se busca optimizar una función objetivo sujeta a ciertas restricciones:\n\n");
    fprintf(f, "\\\\[0.5cm]\n");
    fprintf(f, "\\vspace{0.4cm}\n");
    fprintf(f, "\\textbf{%s} ", opt_type == 0 ? "Maximizar" : "Minimizar \\\\\n");
    fprintf(f, "\\textbf{Z = } ");

    // Nombres de variables (para impresión)
    const char *var_names[64]; // límite práctico: 15; 64 es amplio
    collect_var_names(vars, var_names);

    // Función objetivo
    int first = 1;
    for (int i = 0; i < vars; i++)
    {
        const char *coef = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_obj_coefs, i)));
        const char *var = var_names[i];
        print_term(f, coef, var, first);
        if (fabs(atof(coef)) > 1e-12)
            first = 0;
    }
    fprintf(f, "\\\\[0.5cm]\n");

    // Restricciones
    fprintf(f, "\\textbf{Sujeto a:}\\\\\n");
    int idx = 0;
    for (int r = 0; r < cons; r++)
    {
        int start = 1;
        for (int v = 0; v < vars; v++)
        {
            const char *coef = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_con_coefs, idx++)));
            const char *var = var_names[v];
            print_term(f, coef, var, start);
            if (fabs(atof(coef)) > 1e-12)
                start = 0;
        }
        GtkWidget *combo = g_ptr_array_index(combo_signs, r);
        int sign = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
        const char *sign_text = (sign == 0) ? "\\leq" : (sign == 1) ? "="
                                                                    : "\\geq";
        GtkWidget *rhs = g_ptr_array_index(entry_rhs_list, r);
        const char *rhs_val = gtk_entry_get_text(GTK_ENTRY(rhs));
        fprintf(f, " %s %s \\\\\n", sign_text, rhs_val);
    }

    // No negatividad
    fprintf(f, "\\\\[0.3cm]\n\\textbf{Con } ");
    for (int i2 = 0; i2 < vars; i2++)
    {
        fprintf(f, "%s%s", var_names[i2], (i2 < vars - 1) ? ", " : "");
    }
    fprintf(f, " \\geq 0\\\\[0.5cm]\n\\bigskip\n");

    if (unsupported_signs)
    {
        fprintf(f, "\\textit{Aviso:} En esta versión se resuelven únicamente restricciones de tipo $\\leq$. "
                   "Si se ingresaron restricciones con $=$ o $\\geq$, se ignoran para el proceso de solución y se reportan aquí por completitud.\\\\[0.2cm]\n");
    }

    // --- Tabla inicial ---
    fprintf(f, "\\section*{Tabla inicial}\n");
    latex_print_table(f, run->initial, vars, cons, var_names, -1, -1, 0, -1);

    // --- Tablas intermedias ---
    if (show_tables && run->steps && run->steps->len > 0)
    {
        fprintf(f, "\\section*{Tablas intermedias}\n");
        for (guint k = 0; k < run->steps->len; k++)
        {
            IterStep *step = (IterStep *)g_ptr_array_index(run->steps, k);
            fprintf(f, "\\subsection*{Iteración %u}\n", (unsigned)(k + 1));
            if (step->pre)
            {
                // Mostrar tabla PRE (candidatos y razones)
                int show_ratio = (step->pre->entering_col >= 0);
                latex_print_table(f, step->pre, vars, cons, var_names,
                                  step->pre->entering_col, step->pre->leaving_row,
                                  show_ratio, step->pre->entering_col);
                fprintf(f, "\\vspace{0.2cm}\n");

                // Si hubo degeneración en esta iteración, documentarla (filas empatadas / razón 0)
                if (step->degenerate && step->deg_rows && step->deg_rows->len > 0)
                {
                    fprintf(f, "\\textbf{Degeneración detectada:} razón mínima = %.6g. Filas empatadas o con razón cero:\\\\\n", step->min_ratio);
                    fprintf(f, "\\begin{tabular}{r r r r}\\toprule Fila & Divisor & Dividendo & Razón \\\\ \\midrule\n");
                    for (guint p = 0; p < step->deg_rows->len; p++)
                    {
                        int row = GPOINTER_TO_INT(g_ptr_array_index(step->deg_rows, p));
                        double a = step->pre->T[row][step->pre->entering_col];
                        double rhs = step->pre->T[row][step->pre->cols - 1];
                        double r = (a > EPS) ? rhs / a : 0.0;
                        fprintf(f, "%d & %.6g & %.6g & %.6g \\\\\n", row, a, rhs, r);
                    }
                    fprintf(f, "\\bottomrule\\end{tabular}\\\\[0.2cm]\n");
                    fprintf(f, "Regla usada para romper el empate: se seleccionó la primera fila listada.\\\\[0.2cm]\n");
                }
            }
            if (step->post)
            {
                // Mostrar tabla POST (resultado tras canonizar)
                fprintf(f, "\\textbf{Tabla tras canonizar:}\\\\[0.2cm]\n");
                latex_print_table(f, step->post, vars, cons, var_names,
                                  step->post->entering_col, step->post->leaving_row,
                                  0, -1);
                fprintf(f, "\\vspace{0.3cm}\n");
            }
        }
    }

    // --- Tabla final ---
    fprintf(f, "\\section*{Tabla final}\n");

    if (run->final)
    {
        latex_print_table(f, run->final, vars, cons, var_names,
                          -1, -1, 0, -1);
    }

    if (run->final_alt)
    {
        fprintf(f, "En la solución 1, se puede observar que existe al menos una variable no básica con costo reducido cero, "
                   "lo que indica la presencia de soluciones múltiples óptimas. "
                   "Al realizar un pivote adicional utilizando dicha variable como entrante, se obtiene otra solución básica factible equivalente:\\\\[0.3cm]\n");
        fprintf(f, "\\subsection*{Tabla final 2 (solución 2 tras un pivote adicional)}\n");
        latex_print_table(f, run->final_alt, vars, cons, var_names,
                          -1, -1, 0, -1);
    }

    // --- Soluciones múltiples: imprimir soluciones y combinaciones convexas ---
    fprintf(f, "\\section*{Solución}\n");
    if (run->status == 1)
    {
        /* Explicación extendida para caso no acotado:
           informamos la columna candidata que quedó sin pivotes positivos
           y posibles causas (falta de restricciones o restriciones mal planteadas). */
        fprintf(f, "\\textbf{El problema es no acotado.}\\\\\n");

        /* intentar recuperar la columna entrante reportada en la tabla final */
        int ent_col = -1;
        if (run->final)
            ent_col = run->final->entering_col;

        if (ent_col >= 1)
        {
            if (ent_col <= vars)
            {
                fprintf(f, "La columna candidata para entrar era la columna %d (variable %s).\\\\\n",
                        ent_col, var_names[ent_col - 1]);
            }
            else if (run->final && ent_col < run->final->cols - 1)
            {
                fprintf(f, "La columna candidata para entrar era la columna %d (una columna de holgura).\\\\\n",
                        ent_col);
            }
            else
            {
                fprintf(f, "La columna candidata para entrar era la columna %d.\\\\\n", ent_col);
            }

            fprintf(f,
                    "En esa columna ya no existen coeficientes positivos en las filas de las restricciones (todas las entradas relevantes son \\(\\leq 0\\)), "
                    "por lo que no existe una fila pivote válida que permita continuar el algoritmo.\\\\\n");

            fprintf(f,
                    "Interpretación: la función objetivo puede aumentar (o disminuir, según el sentido) indefinidamente en la dirección asociada a esa variable, "
                    "lo que significa que la región factible no acota el crecimiento de la función objetivo.\\\\\n");

            fprintf(f,
                    "Posibles causas: falta de restricciones que limiten esa dirección (el modelo es incompleto), o restricciones mal planteadas "
                    "(signos o coeficientes incorrectos en las ecuaciones). Revise las restricciones y los signos para corregir la formulación.\\\\[0.2cm]\n");
        }
        else
        {
            fprintf(f, "No se pudo identificar la columna entrante en la tabla final. Revise la entrada de datos y la tabla intermedia para más detalles.\\\\[0.2cm]\n");
        }
    }
    else
    {
        // Solución 1
        /*
        if (run->solution)
        {
            fprintf(f, "$Z = %.6g$\\\\\n", run->z_value);
            fprintf(f, "\\\\\\\n");
            fprintf(f, "\\vspace{1cm}\n");
            fprintf(f, "\\begin{tabular}{lr}\\toprule Variable & Valor \\\\ \\midrule\n");
            for (int j = 0; j < vars; j++)
                fprintf(f, "%s & %.6g \\\\\n", var_names[j], run->solution[j]);
            fprintf(f, "\\bottomrule\\end{tabular}\\\\[0.3cm]\n");
        }
        */

        if (run->final)
        {
            fprintf(f, "$Z = %.6g$\\\\[0.3cm]\n", run->z_value);
            latex_print_all_variables(f, run->final, var_names, vars, cons);
        }
        // Solución 2 (si existe)
        /*
        if (run->solution_alt)
        {

            fprintf(f, "\\subsection*{Solución 2 (tras pivote adicional)}\n");
            fprintf(f, "\\vspace{0.3cm}\n");
            fprintf(f, "\\begin{tabular}{lr}\\toprule Variable & Valor \\\\ \\midrule\n");
            for (int j = 0; j < vars; j++)
                fprintf(f, "%s & %.6g \\\\\n", var_names[j], run->solution_alt[j]);
            fprintf(f, "\\bottomrule\\end{tabular}\\\\[0.3cm]\n");
        }
        */
        if (run->final_alt)
        {
            fprintf(f, "\\subsection*{Solución 2 (tras pivote adicional)}\\\\[0.3cm]\n");
            fprintf(f, "$Z = %.6g$\\\\[0.3cm]\n", run->z_value);
            latex_print_all_variables(f, run->final_alt, var_names, vars, cons);
        }

        // Combinaciones convexas
        if (run->conv_solutions && run->conv_count > 0)
        {
            fprintf(f, "\\subsection*{Soluciones adicionales (combinaciones convexas)}\n");
            fprintf(f, "Si existen soluciones múltiples, cualquier combinación convexa de dos soluciones básicas produce otra solución factible. "
                       "Aquí usamos la fórmula: $x(\\alpha)=\\alpha x^{(1)} + (1-\\alpha) x^{(2)}$, con $\\alpha \\in [0,1]$.\\\\[0.2cm]\n");

            const double alphas_print[3] = {0.25, 0.5, 0.75};
            for (int a = 0; a < run->conv_count; a++)
            {
                fprintf(f, "\\textbf{Alfa = %.2g}:\\\\\n", alphas_print[a]);
                fprintf(f, "\\begin{tabular}{lr}\\toprule Variable & Valor \\\\ \\midrule\n");
                for (int j = 0; j < vars; j++)
                {
                    fprintf(f, "%s & %.6g \\\\\n", var_names[j], run->conv_solutions[a][j]);
                }
                fprintf(f, "\\bottomrule\\end{tabular}\\\\[0.2cm]\n");
            }

            // Mostrar la ecuación general
            fprintf(f, "\\vspace{0.3cm}\n");
        }
    }
    // === GRAFICO PARA CASO DE 2 VARIABLES ===
    if (vars == 2)
    {
        fprintf(f, "\\section*{Gráfico de la región factible (solo para 2 variables)}\n");
        fprintf(f, "A continuación se muestra la región factible en el plano $(x_1, x_2)$ "
                   "junto con las rectas de las restricciones.\\\\[0.4cm]\n");

        // --- PASO 1: Extraer coeficientes ---
        double A[cons][2], B[cons];
        int idx2 = 0;
        for (int r = 0; r < cons; r++)
        {
            A[r][0] = atof(gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_con_coefs, idx2++))));
            A[r][1] = atof(gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_con_coefs, idx2++))));
            B[r] = atof(gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_rhs_list, r))));
        }

        double pts_x[256], pts_y[256];
        int pt_count = 0;

        // Añadir (0,0) siempre como posible punto
        pts_x[pt_count] = 0;
        pts_y[pt_count] = 0;
        pt_count++;

        // Intersección con ejes
        for (int r = 0; r < cons; r++)
        {
            if (fabs(A[r][0]) > 1e-12)
            {
                double x = B[r] / A[r][0];
                if (x >= 0)
                {
                    pts_x[pt_count] = x;
                    pts_y[pt_count] = 0;
                    pt_count++;
                }
            }
            if (fabs(A[r][1]) > 1e-12)
            {
                double y = B[r] / A[r][1];
                if (y >= 0)
                {
                    pts_x[pt_count] = 0;
                    pts_y[pt_count] = y;
                    pt_count++;
                }
            }
        }

        // Intersecciones entre restricciones
        for (int i = 0; i < cons; i++)
            for (int j = i + 1; j < cons; j++)
            {
                double det = A[i][0] * A[j][1] - A[j][0] * A[i][1];
                if (fabs(det) < 1e-12)
                    continue;

                double x = (B[i] * A[j][1] - B[j] * A[i][1]) / det;
                double y = (A[i][0] * B[j] - A[j][0] * B[i]) / det;

                if (x >= 0 && y >= 0)
                {
                    pts_x[pt_count] = x;
                    pts_y[pt_count] = y;
                    pt_count++;
                }
            }

        // Filtrar solo puntos factibles
        double fx[256], fy[256];
        int fc = 0;
        for (int k = 0; k < pt_count; k++)
        {
            int ok = 1;
            for (int r = 0; r < cons; r++)
                if (A[r][0] * pts_x[k] + A[r][1] * pts_y[k] > B[r] + 1e-9)
                    ok = 0;
            if (ok)
            {
                fx[fc] = pts_x[k];
                fy[fc] = pts_y[k];
                fc++;
            }
        }

        // Ordenar vértices por ángulo alrededor del centroide
        double cx = 0, cy = 0;
        for (int i = 0; i < fc; i++)
        {
            cx += fx[i];
            cy += fy[i];
        }
        cx /= fc;
        cy /= fc;

        for (int i = 0; i < fc; i++)
            for (int j = i + 1; j < fc; j++)
            {
                double ang_i = atan2(fy[i] - cy, fx[i] - cx);
                double ang_j = atan2(fy[j] - cy, fx[j] - cx);
                if (ang_i > ang_j)
                {
                    double tx = fx[i], ty = fy[i];
                    fx[i] = fx[j];
                    fy[i] = fy[j];
                    fx[j] = tx;
                    fy[j] = ty;
                }
            }

        // --- DIBUJAR ---
        fprintf(f, "\\begin{tikzpicture}\n");
        fprintf(f, "\\begin{axis}[\n  xlabel=$x_1$, ylabel=$x_2$,\n"
                   "  axis lines=middle,\n  xmin=0, ymin=0,\n"
                   "  width=12cm, height=10cm]\n");

        // Dibujar restricciones completas (rectas enteras)
        for (int r = 0; r < cons; r++)
        {
            if (fabs(A[r][1]) > 1e-12)
                fprintf(f, "\\addplot[thick,domain=0:20] (x,{(%g-%g*x)/%g});\n", B[r], A[r][0], A[r][1]);
            else
            {
                double xval = B[r] / A[r][0];
                fprintf(f, "\\addplot[thick] coordinates {(%g,0) (%g,20)};\n", xval, xval);
            }
        }

        // Región factible (polígono cerrado)
        fprintf(f, "\\addplot[fill=blue!25, fill opacity=0.4] coordinates {");
        for (int i = 0; i < fc; i++)
            fprintf(f, "(%g,%g) ", fx[i], fy[i]);
        fprintf(f, "(%g,%g)};\n", fx[0], fy[0]); // cerrar polígono

        fprintf(f, "\\end{axis}\n\\end{tikzpicture}\n\\newpage\n");
    }

    // --- Información sobre degeneración (si aplica) ---
    if (run->degenerate)
    {
        fprintf(f, "\\subsection*{Problema degenerado}\n");
        fprintf(f,
                "Se detectó degeneración durante la ejecución del Simplex. \\newline\n"
                "Definición: una solución básica factible es degenerada cuando alguna variable básica tiene valor cero. "
                "Esto puede ocurrir cuando en el test de razón mínima la razón mínima es 0 o cuando hay un empate entre razones.\\\\[0.2cm]\n");

        // Explicar cómo se resolvieron empates/degeneración
        fprintf(f,
                "El programa marca como degeneración cuando dos razones difieren en menos de "
                "un epsilon (%.0e) o cuando la razón mínima es 0. Para romper el empate se adopta una heurística simple: se selecciona la primera fila con "
                "la razón mínima encontrada. Esta elección se documenta en las tablas intermedias (fila pivote seleccionada).\\\\[0.2cm]\n",
                EPS);
    }

    // --- Referencias ---
    fprintf(f, "\\renewcommand{\\refname}{Referencias}\n");
    fprintf(f, "\\begin{thebibliography}{9}\n");
    fprintf(f,
            "\\bibitem{wikipedia2025} Wikipedia contributors. (2025, 6 octubre). Simplex algorithm. "
            "\\textit{Wikipedia}. Disponible en: https://en.wikipedia.org/wiki/Simplex_algorithm\n\n");
    fprintf(f,
            "\\bibitem{benlowery2022} Ben-Lowery. (2022, 4 abril). Linear Programming and the birth of the Simplex Algorithm. "
            "\\textit{Ben Lowery @ STOR-i}. Disponible en: https://www.lancaster.ac.uk/stor-i-student-sites/ben-lowery/2022/03/linear-programming-and-the-birth-of-the-simplex-algorithm/\n\n");
    fprintf(f,
            "\\bibitem{libretexts2022} Libretexts. (2022, 18 julio). 4.2: Maximization by the Simplex method. "
            "\\textit{Mathematics LibreTexts}. Disponible en: https://math.libretexts.org/Bookshelves/Applied_Mathematics/Applied_Finite_Mathematics_%%28Sekhon_and_Bloom%%29/04%%3A_Linear_Programming_The_Simplex_Method/4.02%%3A_Maximization_By_The_Simplex_Method\n\n");
    fprintf(f, "\\end{thebibliography}\n");

    fprintf(f, "\\end{document}\n");
    fclose(f);

    tableau_free(initial);
    simplex_run_free(run);
}

// Ejecutar y compilar
void simplex(GtkWidget *btn, gpointer data)
{
    const char *problem_name = g_object_get_data(G_OBJECT(btn), "problem_name");
    if (!problem_name)
        problem_name = "SinNombre";
    GtkWidget **widgets = (GtkWidget **)data;
    GtkWidget *box_model_local = widgets[0];
    int vars = GPOINTER_TO_INT(widgets[1]);
    int cons = GPOINTER_TO_INT(widgets[2]);

    gboolean show_tables = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_intermediate));
    write_file(vars, cons, box_model_local, show_tables, problem_name);

    // Compilar con pdflatex (dos pasadas por si acaso)
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "pdflatex -interaction=nonstopmode %s > /dev/null 2>&1", "simplex_output.tex");
    system(cmd);
    system(cmd);

    // Abrir PDF en evince (presentación)
    snprintf(cmd, sizeof(cmd), "evince --presentation simplex_output.pdf &");
    system(cmd);
}

// --- Main ---
int main(int argc, char *argv[])
{
    gtk_init(&argc, &argv);

    GtkBuilder *builder = gtk_builder_new_from_file("simplex.glade");
    main_window = GTK_WIDGET(gtk_builder_get_object(builder, "main_window"));
    spin_vars = GTK_WIDGET(gtk_builder_get_object(builder, "spin_variables"));
    spin_constraints = GTK_WIDGET(gtk_builder_get_object(builder, "spin_constraints"));
    box_model = GTK_WIDGET(gtk_builder_get_object(builder, "box_model"));

    GtkWidget *btn_generate = GTK_WIDGET(gtk_builder_get_object(builder, "btn_generate"));
    GtkWidget *btn_load = GTK_WIDGET(gtk_builder_get_object(builder, "btn_load"));

    g_signal_connect(btn_load, "clicked", G_CALLBACK(on_load_clicked), builder);
    g_signal_connect(main_window, "destroy", G_CALLBACK(on_close_main_window), NULL);
    g_signal_connect(btn_generate, "clicked", G_CALLBACK(on_generate_clicked), builder);

    gtk_widget_show_all(main_window);
    gtk_main();
    return 0;
}
