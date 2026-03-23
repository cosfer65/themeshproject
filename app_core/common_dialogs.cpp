#include <windows.h>
#include <commdlg.h>

#include "common_dialogs.h"

/**
 * @brief Displays a standard Windows Open File dialog.
 *
 * This function opens a modal file selection dialog, allowing the user to select a file.
 * The dialog is configured to only allow selection of existing files and valid paths.
 *
 * @param file_list A filter string that specifies the file types the dialog displays.
 *                  The format is pairs of display string and filter pattern, separated by '\0', ending with '\0\0'.
 *                  Example: "Text Files\0*.TXT\0All Files\0*.*\0"
 * @return const char* Returns a pointer to a static buffer containing the selected file path if successful,
 *                     or nullptr if the dialog is canceled or an error occurs.
 *
 * @note The returned pointer is valid until the next call to this function.
 * @note This function is not thread-safe due to the use of a static buffer.
 */
const char* OpenFileDialog(const char* file_list) {
    OPENFILENAME ofn; // Common dialog box structure
    static char szFile[1024] = { 0 }; // Buffer for file name
    //memset(szFile, 0, sizeof(szFile));
    szFile[0] = 0;
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.hwndOwner = NULL; // Handle to owner window
    ofn.lpstrFile = szFile;
    ofn.nMaxFile = sizeof(szFile);
    ofn.lpstrFilter = file_list;
    ofn.nFilterIndex = 1;
    ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
    if (GetOpenFileName(&ofn)) {
        return szFile;
    }
    return nullptr;
}
