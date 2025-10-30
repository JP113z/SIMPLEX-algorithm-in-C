#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// ========================
//      CONSTANTES
// ========================
#define EPS 1e-9
#define MAX_ITERS 10000

// ========================
//  WIDGETS (tuyos)
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

// --- Prototipos GUI (tuyos) ---
void on_generate_clicked(GtkButton *button, gpointer user_data);
void simplex(GtkWidget *btn, gpointer data);
void write_file(int vars, int cons, GtkWidget *box_model, int show_tables, const char *problem_name);

// ========================
// === NUEVO: Estructuras
// ========================
typedef struct {
    int rows;      // 1 (fila objetivo) + cons
    int cols;      // 1 (Z) + vars + cons (holguras) + 1 (RHS)
    double **T;    // matriz [rows][cols]
    int *basis;    // índice de columna básica por cada fila (excepto fila 0). basis[r] = columna básica de la fila r (r>=1)
    int entering_col; // última col que entró
    int leaving_row;  // última fila que salió
} Tableau;

typedef struct {
    GPtrArray *steps; // Tableau* de cada iteración (se guardan copias)
    Tableau *initial; // copia de la tabla inicial
    Tableau *final;   // copia de la tabla final
    int status;       // 0=Óptimo, 1=No acotado, 2=No factible (no usamos en esta versión), 3=Degenerado (bandera informativa)
    int has_multiple; // 1 si hay soluciones múltiples
    int degenerate;   // 1 si hubo razón mínima 0 o empates (informativo)
    double *solution; // valores de variables de decisión (longitud = vars)
    double z_value;   // valor óptimo de Z
} SimplexRun;

// ========================
// === NUEVO: Utilidades
// ========================
static double **alloc_matrix(int r, int c) {
    double **m = (double**)malloc(sizeof(double*) * r);
    for (int i = 0; i < r; i++) {
        m[i] = (double*)calloc(c, sizeof(double));
    }
    return m;
}
static void free_matrix(double **m, int r) {
    if (!m) return;
    for (int i = 0; i < r; i++) free(m[i]);
    free(m);
}

static Tableau* tableau_new(int rows, int cols) {
    Tableau *tb = (Tableau*)calloc(1, sizeof(Tableau));
    tb->rows = rows;
    tb->cols = cols;
    tb->T = alloc_matrix(rows, cols);
    tb->basis = (int*)malloc(sizeof(int) * rows);
    for (int i = 0; i < rows; i++) tb->basis[i] = -1;
    tb->entering_col = -1;
    tb->leaving_row  = -1;
    return tb;
}
static Tableau* tableau_copy(const Tableau* src) {
    Tableau *t = tableau_new(src->rows, src->cols);
    for (int i = 0; i < src->rows; i++)
        memcpy(t->T[i], src->T[i], sizeof(double)*src->cols);
    memcpy(t->basis, src->basis, sizeof(int)*src->rows);
    t->entering_col = src->entering_col;
    t->leaving_row = src->leaving_row;
    return t;
}
static void tableau_free(Tableau *tb) {
    if (!tb) return;
    free_matrix(tb->T, tb->rows);
    free(tb->basis);
    free(tb);
}

static SimplexRun* simplex_run_new(void) {
    SimplexRun *r = (SimplexRun*)calloc(1, sizeof(SimplexRun));
    r->steps = g_ptr_array_new_with_free_func((GDestroyNotify)tableau_free);
    r->status = 0;
    r->has_multiple = 0;
    r->degenerate = 0;
    r->solution = NULL;
    r->z_value = 0.0;
    return r;
}
static void simplex_run_free(SimplexRun *run) {
    if (!run) return;
    if (run->initial) tableau_free(run->initial);
    if (run->final) tableau_free(run->final);
    if (run->steps) g_ptr_array_free(run->steps, TRUE);
    if (run->solution) free(run->solution);
    free(run);
}

// ========================
// === NUEVO: Construcción de la tabla inicial
// Forma: [fila 0] Z - Σ c_i x_i - Σ 0*s_i = 0
// Restricciones ≤: Σ a_ij x_j + s_i = b_i
// Columnas: [Z | x_1..x_n | s_1..s_m | RHS]
// Filas:    [0 |     1..n |  1..m    |  0 ] (fila 0 con -c en x_j)
// ========================
static Tableau* build_initial_tableau_from_inputs(int vars, int cons, int *unsupported_sign_found, int opt_type) {
    int rows = cons + 1;
    int cols = 1 /*Z*/ + vars + cons /*slacks*/ + 1 /*RHS*/;

    Tableau *tb = tableau_new(rows, cols);

    // Fila 0: Z - sum c x = 0
    tb->T[0][0] = 1.0; // coeficiente de Z
    for (int j = 0; j < vars; j++) {
        const char *coef_str = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_obj_coefs, j)));
        double c = atof(coef_str);
        // En la fila 0 ponemos -c (forma estándar)
        tb->T[0][1 + j] = -c;
    }
    // RHS en fila 0 = 0

    // Restricciones
    *unsupported_sign_found = 0;
    int idx = 0;
    for (int i = 0; i < cons; i++) {
        // Coeficientes x_j
        for (int j = 0; j < vars; j++) {
            const char *a_str = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_con_coefs, idx++)));
            double a = atof(a_str);
            tb->T[i+1][1 + j] = a;
        }
        // Signo
        GtkWidget *combo = g_ptr_array_index(combo_signs, i);
        int sign = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
        if (sign != 0) {
            // En esta versión, solo ≤. Marcamos que hay restricciones no soportadas
            *unsupported_sign_found = 1;
        }
        // Holgura s_i
        tb->T[i+1][1 + vars + i] = 1.0;

        // RHS
        GtkWidget *rhs = g_ptr_array_index(entry_rhs_list, i);
        const char *rhs_str = gtk_entry_get_text(GTK_ENTRY(rhs));
        double b = atof(rhs_str);
        tb->T[i+1][cols - 1] = b;

        // Base inicial: s_i
        tb->basis[i+1] = 1 + vars + i; // columna de s_i
    }

    // Para MINIMIZACIÓN según indicación del usuario:
    // "En casos de minimización se elige la columna a canonizar del valor más positivo de la primera fila"
    // -> mantenemos la fila 0 como -c. La regla de parada/selección se ajusta en la iteración, no aquí.
    (void)opt_type;

    return tb;
}

// ========================
// === NUEVO: Selección de columna que entra
//   max: más NEGATIVO en fila 0 (entre x y s no básicos, normalmente x)
//   min: más POSITIVO en fila 0
// Retorna -1 si ya óptimo según el criterio
// ========================
static int choose_entering_col(const Tableau *tb, int vars, int cons, int is_max) {
    int start = 1; // después de Z
    int end   = 1 + vars + cons - 1; // última col de slack; RHS es la última col
    double best_val = is_max ? 1e100 : -1e100; // max: buscamos el más negativo => inicial alto; min: más positivo => inicial bajo
    int best_col = -1;

    for (int j = start; j < tb->cols - 1; j++) {
        double v = tb->T[0][j];
        if (is_max) {
            if (v < best_val - EPS) { // más negativo
                best_val = v;
                best_col = j;
            }
        } else { // min
            if (v > best_val + EPS) { // más positivo
                best_val = v;
                best_col = j;
            }
        }
    }

    // Criterio de optimalidad:
    // max: si NO hay negativos => óptimo
    // min: si NO hay positivos => óptimo
    if (is_max) {
        if (best_col == -1 || best_val >= -EPS) return -1;
    } else {
        if (best_col == -1 || best_val <= EPS) return -1;
    }
    return best_col;
}

// ========================
// === NUEVO: Selección de fila que sale (test de razón mínima)
// Solo razones positivas. Marca degeneración si hay 0 o empates cercanos
// ========================
static int choose_leaving_row(const Tableau *tb, int entering_col, int *degenerate_flag) {
    int best_row = -1;
    double best_ratio = 0.0;
    *degenerate_flag = 0;

    for (int i = 1; i < tb->rows; i++) {
        double a = tb->T[i][entering_col];
        if (a > EPS) { // solo positivos
            double rhs = tb->T[i][tb->cols - 1];
            double ratio = rhs / a;
            if (best_row == -1 || ratio < best_ratio - EPS) {
                best_ratio = ratio;
                best_row = i;
            } else if (fabs(ratio - best_ratio) <= EPS) {
                *degenerate_flag = 1; // empate de fracciones
            }
            if (ratio <= EPS) *degenerate_flag = 1; // degeneración (sale con 0)
        }
    }
    return best_row; // -1 => No acotado
}

// ========================
// === NUEVO: Pivoteo (Gauss-Jordan en la tabla)
// ========================
static void pivot(Tableau *tb, int pr, int pc) {
    double piv = tb->T[pr][pc];
    // Normalizar fila pivote
    for (int j = 0; j < tb->cols; j++) tb->T[pr][j] /= piv;

    // Anular columna pivote en las demás filas
    for (int i = 0; i < tb->rows; i++) {
        if (i == pr) continue;
        double factor = tb->T[i][pc];
        if (fabs(factor) > EPS) {
            for (int j = 0; j < tb->cols; j++) {
                tb->T[i][j] -= factor * tb->T[pr][j];
            }
        }
    }
    tb->basis[pr] = pc;
    tb->entering_col = pc;
    tb->leaving_row  = pr;
}

// ========================
// === NUEVO: Detección de soluciones múltiples
// max: en óptimo, si existe costo reducido 0 en columna NO básica (x), hay múltiples
// min: análogo (fila 0 con 0 en col no básica)
// ========================
static int detect_multiple_solutions(const Tableau *tb) {
    // revisamos columnas no básicas en fila 0
    for (int j = 1; j < tb->cols - 1; j++) {
        // si la columna no está en la base:
        int in_basis = 0;
        for (int i = 1; i < tb->rows; i++) if (tb->basis[i] == j) { in_basis = 1; break; }
        if (!in_basis) {
            if (fabs(tb->T[0][j]) <= EPS) return 1;
        }
    }
    return 0;
}

// ========================
// === NUEVO: Ejecutar Símplex (max/min)
// Guarda pasos si show_steps != 0
// ========================
static SimplexRun* simplex_solve(Tableau *init, int vars, int cons, int is_max, int show_steps) {
    SimplexRun *run = simplex_run_new();
    run->initial = tableau_copy(init);
    Tableau *tb = tableau_copy(init);

    if (show_steps) {
        g_ptr_array_add(run->steps, tableau_copy(tb));
    }

    int iters = 0;
    while (iters++ < MAX_ITERS) {
        int entering = choose_entering_col(tb, vars, cons, is_max);
        if (entering == -1) { // óptimo
            run->status = 0;
            break;
        }
        int deg_flag = 0;
        int leaving = choose_leaving_row(tb, entering, &deg_flag);
        if (leaving == -1) { // no acotado
            run->status = 1;
            tb->entering_col = entering;
            tb->leaving_row  = -1;
            if (show_steps) g_ptr_array_add(run->steps, tableau_copy(tb));
            break;
        }
        if (deg_flag) run->degenerate = 1;

        pivot(tb, leaving, entering);
        if (show_steps) g_ptr_array_add(run->steps, tableau_copy(tb));
    }

    run->final = tableau_copy(tb);

    // solución y Z
    run->solution = (double*)calloc(vars, sizeof(double));
    for (int i = 1; i < tb->rows; i++) {
        int bc = tb->basis[i];
        // si bc corresponde a alguna x_j
        if (bc >= 1 && bc <= vars) {
            run->solution[bc - 1] = tb->T[i][tb->cols - 1];
        }
    }
    // Valor de Z:
    // Nuestra forma es Z - sum c x = 0 -> RHS fila 0 = Z
    run->z_value = tb->T[0][tb->cols - 1];

    // múltiples soluciones si óptimo y hay costo reducido 0 en alguna no básica
    if (run->status == 0) {
        run->has_multiple = detect_multiple_solutions(tb);
        // Si hay múltiples, "pivoteamos una vez más" a otra solución básica alterna
        // (simple heurística: escoger una col no básica con costo reducido ~ 0 y tratar de pivotear si posible)
        if (run->has_multiple) {
            for (int j = 1; j < tb->cols - 1; j++) {
                int in_basis = 0;
                for (int i = 1; i < tb->rows; i++) if (tb->basis[i] == j) { in_basis = 1; break; }
                if (!in_basis && fabs(tb->T[0][j]) <= EPS) {
                    int dummy_deg = 0;
                    int lv = choose_leaving_row(tb, j, &dummy_deg);
                    if (lv != -1) {
                        pivot(tb, lv, j);
                        break;
                    }
                }
            }
        }
    }

    tableau_free(tb);
    return run;
}

// ========================
// === NUEVO: Impresores LaTeX
// ========================
static void latex_print_table(FILE *f, const Tableau *tb, int vars, int cons,
                              const char **var_names, int highlight_col, int highlight_row,
                              int show_ratio, int ratio_col) {
    // Encabezados: Z | x's | s's | RHS
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "\\rowcolors{2}{white}{gray!10}\n");
    fprintf(f, "\\begin{tabular}{");
    for (int j = 0; j < tb->cols; j++) fprintf(f, "r");
    fprintf(f, "}\n\\toprule\n");

    // Primera fila de nombres
    fprintf(f, "$Z$");
    for (int j = 0; j < vars; j++) fprintf(f, " & %s", var_names[j]);
    for (int j = 0; j < cons; j++) fprintf(f, " & $s_{%d}$", j+1);
    fprintf(f, " & RHS\\\\\\midrule\n");

    // Filas
    for (int i = 0; i < tb->rows; i++) {
        for (int j = 0; j < tb->cols; j++) {
            int is_pivot = (i == highlight_row && j == highlight_col);
            if (is_pivot) fprintf(f, "\\cellcolor{yellow!40}");
            fprintf(f, "%s%.6g", (j==0?"": " & "), tb->T[i][j]);
        }
        fprintf(f, "\\\\\n");
    }
    fprintf(f, "\\bottomrule\n\\end{tabular}\n");

    // Mostrar razones si se pide (columna ratio_col)
    if (show_ratio && ratio_col >= 0) {
        fprintf(f, "\\\\[0.2cm]\n\\small Razones $\\frac{RHS}{a_{i,%d}}$ (sólo con $a_{i,%d}>0$):\\\\\n", ratio_col, ratio_col);
        fprintf(f, "\\begin{tabular}{lr}\\toprule Fila & Razón \\\\ \\midrule\n");
        for (int i = 1; i < tb->rows; i++) {
            double a = tb->T[i][ratio_col];
            if (a > EPS) {
                double rhs = tb->T[i][tb->cols - 1];
                double r = rhs / a;
                fprintf(f, "%d & %.6g \\\\\n", i, r);
            }
        }
        fprintf(f, "\\bottomrule\\end{tabular}\n");
    }

    fprintf(f, "\\end{center}\n");
}

static void latex_print_solution(FILE *f, const SimplexRun *run, int vars, const char **var_names) {
    fprintf(f, "\\subsection*{Solución óptima}\n");
    if (run->status == 1) {
        fprintf(f, "\\textbf{El problema es no acotado.}\\\\\n");
        return;
    }
    fprintf(f, "$Z^* = %.6g$\\\\\n", run->z_value);
    fprintf(f, "Valores de las variables:\\\\\n");
    fprintf(f, "\\begin{tabular}{lr}\\toprule Variable & Valor \\\\ \\midrule\n");
    for (int j = 0; j < vars; j++) {
        fprintf(f, "%s & %.6g \\\\\n", var_names[j], run->solution ? run->solution[j] : 0.0);
    }
    fprintf(f, "\\bottomrule\\end{tabular}\\\\[0.3cm]\n");

    if (run->has_multiple) {
        fprintf(f, "\\textbf{Soluciones múltiples:} Se detectó al menos una dirección con costo reducido cero.\\newline\n");
        fprintf(f, "La ecuación paramétrica de soluciones se obtiene dejando libre la variable no básica con costo reducido cero y ajustando la base.\\newline\n");
        fprintf(f, "Se muestran al menos dos soluciones básicas factibles equivalentes (ver tablas finales).\\\\\n");
    }
    if (run->degenerate) {
        fprintf(f, "\\textit{Observación:} Se detectó \\textbf{degeneración} (razones con empate o salidas con RHS=0).\\\\\n");
    }
}

// ========================
//  TUS HANDLERS GUI (sin cambios salvo donde integramos Símplex)
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

    g_signal_connect_data(btn_exec, "clicked", G_CALLBACK(simplex), data,
                          (GClosureNotify)g_free, 0);

    for (int i = 0; i < entry_var_names->len; i++)
    {
        GtkWidget *entry = g_ptr_array_index(entry_var_names, i);
        g_signal_connect(entry, "changed", G_CALLBACK(on_var_name_changed), GINT_TO_POINTER(i));
    }

    gtk_widget_show_all(second_window);
}

// === NUEVO: helper para obtener nombres de variables en arreglo C
static void collect_var_names(int vars, const char **names_buf) {
    for (int i = 0; i < vars; i++) {
        names_buf[i] = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, i)));
    }
}

// --- Escribimos TODO el PDF incluyendo Simplex ---
void write_file(int vars, int cons, GtkWidget *box_model, int show_tables, const char *problem_name)
{
    // 1) Construir tabla inicial del Símplex
    int unsupported_signs = 0;
    int opt_type = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_opt_global)); // 0=max, 1=min
    Tableau *initial = build_initial_tableau_from_inputs(vars, cons, &unsupported_signs, opt_type);

    // 2) Resolver
    SimplexRun *run = simplex_solve(initial, vars, cons, opt_type==0, show_tables!=0);

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
    fprintf(f, "\\geometry{margin=0.8in}\n");
    fprintf(f, "\\title{Proyecto 4: SIMPLEX}\n");
    fprintf(f, "\\author{Investigación de Operaciones}\n");
    fprintf(f, "\\date{}\n");
    fprintf(f, "\\begin{document}\n");

    // --- Portada (dejé la tuya; puedes ajustar nombres/profesor) ---
    fprintf(f, "\\begin{center}\n");
    fprintf(f, "{\\large Instituto Tecnológico de Costa Rica\\\\[1cm]\n");
    fprintf(f, "\\includegraphics[width=0.4\\textwidth]{TEC.png}\\\\[2cm]\n");
    fprintf(f, "{\\LARGE \\textbf{Proyecto 4: SIMPLEX}}\\\\[2cm]\n");
    fprintf(f, "{\\large Investigación de Operaciones\\\\[2cm]\n");
    fprintf(f, "{\\large Profesor: }\\\\[1cm]\n");
    fprintf(f, "{\\large Francisco Jose Torres Roja}\\\\[2cm]\n");
    fprintf(f, "{\\large Integrantes: }\\\\[1cm]\n");
    fprintf(f, "{\\large Jose Pablo Fernandez Jimenez - 2023117752}\\\\[1cm]\n");
    fprintf(f, "{\\large Diego Durán Rodríguez - 2022437509}\\\\[2cm]\n");
    fprintf(f, "{\\large Segundo semestre 2025\\\\[1cm]\n");
    fprintf(f, "\\end{center}\n\\newpage\n");

    // --- Algoritmo Simplex (tu texto) ---
    fprintf(f, "\\section*{Algoritmo SIMPLEX}\n");
    fprintf(f,
            "El \\textbf{algoritmo SIMPLEX} resuelve problemas de programación lineal avanzando por vértices del poliedro factible hasta hallar el óptimo. "
            "En este trabajo seguimos la forma canónica: en la tabla inicial la primera fila corresponde a la función objetivo $Z - \\sum c_i x_i = 0$, "
            "y cada fila siguiente a una restricción con su variable de holgura. En maximización se elige como columna entrante el \\emph{más negativo} de la fila 0; "
            "en minimización, el \\emph{más positivo} (como se indicó en clase). El pivote se selecciona por \\emph{mínima razón positiva} y se documentan los casos especiales.\n");

    // --- Problema original ---
    fprintf(f, "\\section*{Problema original}\n");
    fprintf(f, "Nombre del problema: \\textbf{%s}\\\\[0.3cm]\n", problem_name);
    fprintf(f, "\\textbf{%s Z = } ", opt_type == 0 ? "Maximizar" : "Minimizar");

    // Nombres de variables (para impresión)
    const char *var_names[64]; // límite práctico: 15; 64 es amplio
    collect_var_names(vars, var_names);

    // Función objetivo
    for (int i = 0; i < vars; i++) {
        const char *coef = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_obj_coefs, i)));
        const char *var = var_names[i];
        if (i > 0) fprintf(f, " + ");
        fprintf(f, "%s %s", coef, var);
    }
    fprintf(f, "\\\\[0.5cm]\n");

    // Restricciones
    fprintf(f, "\\textbf{Sujeto a:}\\\\\n");
    int idx = 0;
    for (int r = 0; r < cons; r++) {
        for (int v = 0; v < vars; v++) {
            const char *coef = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_con_coefs, idx++)));
            const char *var = var_names[v];
            if (v > 0) fprintf(f, " + ");
            fprintf(f, "%s %s", coef, var);
        }
        GtkWidget *combo = g_ptr_array_index(combo_signs, r);
        int sign = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
        const char *sign_text = (sign == 0) ? "\\leq" : (sign == 1) ? "=" : "\\geq";
        GtkWidget *rhs = g_ptr_array_index(entry_rhs_list, r);
        const char *rhs_val = gtk_entry_get_text(GTK_ENTRY(rhs));
        fprintf(f, " %s %s \\\\\n", sign_text, rhs_val);
    }

    // No negatividad
    fprintf(f, "\\\\[0.3cm]\n\\textbf{Y con } ");
    for (int i2 = 0; i2 < vars; i2++) {
        fprintf(f, "%s%s", var_names[i2], (i2 < vars - 1) ? ", " : "");
    }
    fprintf(f, " \\geq 0\\\\[0.5cm]\n\\bigskip\n");

    if (unsupported_signs) {
        fprintf(f, "\\textit{Aviso:} En esta versión se resuelven únicamente restricciones de tipo $\\leq$. "
                    "Si se ingresaron restricciones con $=$ o $\\geq$, se ignoran para el proceso de solución y se reportan aquí por completitud.\\\\[0.2cm]\n");
    }

    // --- Tabla inicial ---
    fprintf(f, "\\section*{Tabla inicial}\n");
    latex_print_table(f, run->initial, vars, cons, var_names, -1, -1, 0, -1);

    // --- Tablas intermedias ---
    if (show_tables && run->steps && run->steps->len > 0) {
        fprintf(f, "\\section*{Tablas intermedias}\n");
        for (guint k = 0; k < run->steps->len; k++) {
            Tableau *step = (Tableau*)g_ptr_array_index(run->steps, k);
            fprintf(f, "\\subsection*{Iteración %u}\n", (unsigned)k);
            int show_ratio = (step->entering_col >=0);
            latex_print_table(f, step, vars, cons, var_names,
                              step->entering_col, step->leaving_row,
                              show_ratio, step->entering_col);
            fprintf(f, "\\vspace{0.3cm}\n");
        }
    }

    // --- Tabla final ---
    fprintf(f, "\\section*{Tabla final}\n");
    latex_print_table(f, run->final, vars, cons, var_names,
                      -1, -1, 0, -1);

    // --- Solución ---
    fprintf(f, "\\section*{Solución}\n");
    latex_print_solution(f, run, vars, var_names);

    // --- Referencias (tuyas) ---
    fprintf(f, "\\section*{Referencias}\n");
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
    // Nota: si recibes "archivo no existe", revisa permisos/dir trabajo o errores de pdflatex.
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
