#NoEnv  ; Recommended for performance and compatibility with future AutoHotkey releases.
; #Warn  ; Enable warnings to assist with detecting common errors.
SendMode Input  ; Recommended for new scripts due to its superior speed and reliability.
SetWorkingDir %A_ScriptDir%  ; Ensures a consistent starting directory.

;This script modifies the arrow keys inside (and only inside) Adobe Reader DC to perform scrolling


SetTitleMatchMode 2 ;Makes the WinTitle parameter, used by IfWinActive, select a substring of the title of a window

#IfWinActive, Adobe Acrobat Reader DC
Up::
Send {WheelUp 5}
return


#IfWinActive, Adobe Acrobat Reader DC
Down::
Send {WheelDown 5}
return

#IfWinActive, Adobe Acrobat Reader DC
XButton1::
Send !{Left}
return

#IfWinActive, Adobe Acrobat Reader DC
XButton2::
Send !{Right}
return

WheelUp::
Send {WheelUp 4}
return

WheelDown::
Send {WheelDown 4}
return


