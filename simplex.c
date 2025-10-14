#include <gtk/gtk.h>
#include <stdio.h>
#include <string.h>

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

// --- Prototipos ---
void on_generate_clicked(GtkButton *button, gpointer user_data);
void simplex(GtkWidget *btn, gpointer data);
void write_file(int vars, int cons, GtkWidget *box_model, int show_tables, const char *problem_name);

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

void write_file(int vars, int cons, GtkWidget *box_model, int show_tables, const char *problem_name)
{
    FILE *f = fopen("simplex_output.tex", "w");
    if (!f)
    {
        perror("fopen tex");
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
    fprintf(f, "\\begin{document}\n");

    // --- Portada ---
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

    // --- Algoritmo Simplex ---
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
    fprintf(f, "El problema original se puede formular como un problema de programación lineal, donde se busca optimizar una función objetivo sujeta a ciertas restricciones.\n\n");

    // Obtener tipo de optimización
    int opt_type = gtk_combo_box_get_active(GTK_COMBO_BOX(combo_opt_global));
    fprintf(f, "\\textbf\n\n{%s Z = } ", opt_type == 0 ? "Maximizar" : "Minimizar");

    // --- Escribir función objetivo ---
    for (int i = 0; i < vars; i++)
    {
        const char *coef = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_obj_coefs, i)));
        const char *var = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, i)));

        if (i > 0)
            fprintf(f, " + ");
        fprintf(f, "%s %s", coef, var);
    }

    fprintf(f, "\\\\[0.5cm]\n");

    // --- Escribir restricciones ---
    fprintf(f, "\\textbf{Sujeto a:}\\\\\n");

    int idx = 0;
    for (int r = 0; r < cons; r++)
    {
        for (int v = 0; v < vars; v++)
        {
            const char *coef = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_con_coefs, idx++)));
            const char *var = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, v)));

            if (v > 0)
                fprintf(f, " + ");
            fprintf(f, "%s %s", coef, var);
        }

        GtkWidget *combo = g_ptr_array_index(combo_signs, r);
        int sign = gtk_combo_box_get_active(GTK_COMBO_BOX(combo));
        const char *sign_text = (sign == 0) ? "\\leq" : (sign == 1) ? "="
                                                                    : "\\geq";

        GtkWidget *rhs = g_ptr_array_index(entry_rhs_list, r);
        const char *rhs_val = gtk_entry_get_text(GTK_ENTRY(rhs));

        fprintf(f, " %s %s \\\\\n", sign_text, rhs_val);
    }

    // --- Escribir condición de no negatividad ---
    fprintf(f, "\\\\[0.3cm]\n");
    fprintf(f, "\\textbf{Y con } ");

    for (int i = 0; i < vars; i++)
    {
        const char *var = gtk_entry_get_text(GTK_ENTRY(g_ptr_array_index(entry_var_names, i)));
        fprintf(f, "%s", var);
        if (i < vars - 1)
            fprintf(f, ", ");
    }

    fprintf(f, " \\geq 0\\\\[0.5cm]\n\\bigskip\n");

    // --- Tabla inicial---
    fprintf(f, "\\section*{Tabla inicial}\n");

    // --- Tablas intermedias ---
    if (show_tables)
    {
        fprintf(f, "\\section*{Tablas intermedias}\n");
    }

    fprintf(f, "\\newpage\n");

    // --- Tabla final ---
    fprintf(f, "\\section*{Tabla final}\n");

    // --- Soluciones ---
    fprintf(f, "\\section*{Soluciones}\n");

    // --- Referencias ---
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
}

void simplex(GtkWidget *btn, gpointer data)
{
    const char *problem_name = g_object_get_data(G_OBJECT(btn), "problem_name");
    if (!problem_name)
        problem_name = "SinNombre";
    GtkWidget **widgets = (GtkWidget **)data;
    GtkWidget *box_model = widgets[0];
    int vars = GPOINTER_TO_INT(widgets[1]);
    int cons = GPOINTER_TO_INT(widgets[2]);

    gboolean show_tables = gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(check_intermediate));
    write_file(vars, cons, box_model, show_tables, problem_name);

    // Compilar con pdflatex
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
