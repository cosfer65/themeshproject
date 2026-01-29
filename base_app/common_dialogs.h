#ifndef __common_dialogs__
#define __common_dialogs__

/**
 * @brief Displays an open file dialog to the user.
 *
 * Presents a dialog box that allows the user to select a file from the file system.
 * The dialog filters files according to the specified file list.
 *
 * @param file_list A string specifying the file types to filter (e.g., "*.txt;*.png").
 * @return A pointer to a string containing the selected file path, or nullptr if no file was selected.
 */
const char* OpenFileDialog(const char* file_list);

#endif // __common_dialogs__
