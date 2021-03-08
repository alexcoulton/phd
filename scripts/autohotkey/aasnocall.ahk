#NoEnv  ; Recommended for performance and compatibility with future AutoHotkey releases.
; #Warn  ; Enable warnings to assist with detecting common errors.
SendMode Input  ; Recommended for new scripts due to its superior speed and reliability.
SetWorkingDir %A_ScriptDir%  ; Ensures a consistent starting directory.







Numpad0::
BlockInput On
MouseGetPos, xpos, ypos
Click Right
MouseMove, (xpos+30), (ypos+40)
MouseGetPos, xpos, ypos
Sleep, 700
MouseMove, (xpos+300), (ypos+100)
Sleep, 200
Click Left
Sleep, 500
MouseMove, 500, 1300
Click Left
MouseMove, xpos, ypos
BlockInput Off

return