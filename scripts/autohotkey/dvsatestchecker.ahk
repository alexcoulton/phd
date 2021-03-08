#NoEnv  ; Recommended for performance and compatibility with future AutoHotkey releases.
; #Warn  ; Enable warnings to assist with detecting common errors.
SendMode Input  ; Recommended for new scripts due to its superior speed and reliability.
SetWorkingDir %A_ScriptDir%  ; Ensures a consistent starting directory.

Numpad1::
#IfWinExist, AxC 51-2 - alex
	WinActivate ;

WinGetTitle, Title, A
;MsgBox, The active window is "%Title%".
; MsgBox, blablablablablablablab
MouseGetPos, xpos, ypos	
MsgBox, "%xpos%" "%ypos%"
Run "C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"
Run https://driverpracticaltest.dvsa.gov.uk/login				
Send {tab}
Send {tab}
Send {tab}
Send {tab}

return


Numpad2::
MouseGetPos, xpos, ypos
MsgBox, blablablablaba2
MsgBox, "%xpos%" "%ypos%"
return

Numpad3::
Send {tab}
return


seqone(textto1)
{
	MouseMove, 526, 137
	MouseClick
	Sleep, 100
	MouseMove, 533, 173
	MouseClick
	Sleep, 100
	MouseMove 188, 131
	MouseClick
	Sleep, 100
	MouseMove, 124, 154
	MouseClick
	MouseMove, 370, 130
	MouseClick
	MouseMove, 342, 193
	MouseClick
	MouseMove, 459, 133
	MouseClick
	Send ^a
	Send %Clipboard%
	MouseMove, 464, 474
	MouseClick
	Sleep, 100
	MouseMove, 103, 243
	MouseClick
	return
}

seqtwo()
{
	MouseMove, 74, 1001
	MouseClick
	Send, ^a
	Send, %Clipboard%
	MouseMove, 191, 998
	MouseClick
}

saveimage(num1)
{
	MouseMove, 1260, 94
	MouseClick
	Sleep, 500
	Send %Clipboard%	
	Send %num1%
	Sleep, 500
	MouseMove, 762, 863
	MouseClick
}





; Numpad4::
; #IfWinExist, c.x.a65-3 - alex
; 	WinActivate 
; seqtwo()
; Sleep, 500
; #IfWinExist, c.x.a65-1 - alex
; 	WinActivate 
; seqtwo()
; Sleep, 500
; #IfWinExist, AxC 61-1 - alex
; 	WinActivate 
; seqtwo()
; Sleep, 500
; #IfWinExist, AxC 51-2 - alex
; 	WinActivate 
; seqtwo()
; Sleep, 500
; return

Numpad4::
#IfWinExist, c.x.a65-3 - alex
	WinActivate, c.x.a65-3 - alex
seqtwo()
Sleep, 200
saveimage("cxa653")
Sleep, 500

#IfWinExist, c.x.a65-1 - alex
	WinActivate, c.x.a65-1 - alex
seqtwo()
Sleep, 200
saveimage("cxa651")
Sleep, 500

#IfWinExist, AxC 61-1 - alex
	WinActivate, AxC 61-1 - alex
seqtwo()
Sleep, 200
saveimage("axc611")
Sleep, 500

#IfWinExist, AxC 51-2 - alex
	WinActivate, AxC 51-2 - alex 
seqtwo()
Sleep, 200
saveimage("axc512")
Sleep, 500
return