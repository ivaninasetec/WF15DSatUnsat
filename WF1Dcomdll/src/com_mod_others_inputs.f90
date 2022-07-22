	!********************************************************************************************************************
	! TITLE         : LIBRARY OF FUNCTIONS TO USE IN THE INPUT PROCESS
	! PROJECT       : WF1DCOMDLL
	! MODULE        : com_mod_inputs
	! URL           : https://github.com/ivaninasetec/WF15DSatUnsat
	! AFFILIATION   : The University of Nottingham
	! DATE          : 13/2/2022
	! REVISION      : 1.0
	! LICENSE       : This software is copyrighted 2022(C)
	!
	! DESCRIPTION:
	!> Class for the collection of nodes-classes in the saturated model
	!>
	!> @author
	!> Iván Campos-Guereta Díez
	!> MSc Civil Engineering by <a href="http://www.upm.es/">Polytechnic University of Madrid</a>
	!> PhD Student by <a href="https://www.nottingham.ac.uk/">The university of Nottingham</a>
	!> eMBA by <a href="https://www.santelmo.org/en">San Telmo Bussiness School</a>
	!> ivan.camposguereta@nottingham.ac.uk
	!> Working partner of <a href="https://www.inasetec.es">INASETEC</a>
	!********************************************************************************************************************

	module com_mod_inputs

	implicit none
	include 'inc_precision.fi'

	private

	public::s_com_inputs_locateblock,s_com_inputs_nextrecord
	contains

	!********************************************************************************************************************
	! S: S_COM_INPUTS_LOCATEBLOCK(FileIndex, Bloq, Readerror)
	!--------------------------------------------------------------------------------------------------------------------
	! This subroutine rewind the file with index "FileIndex" and search for the line with the word: 'Block '+Bloq
	! that will define the initial line for the input of that block.
	! In case of error returns Readerror<>0
	!********************************************************************************************************************

	subroutine s_com_inputs_locateblock(fileindex, bloq, readerror)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_inputs_locateblock" :: s_com_inputs_locateblock
	!DEC$ endif

	! search for 'block x' in a file and let the file pointer at that row

	integer, intent(in):: fileindex
	character*1, intent(in)::bloq
	integer, intent(out)::readerror

	character*400:: chblock
	logical::checkblock

	rewind fileindex
	readerror = 0
	checkblock=.false.
	do while (.not.checkblock)
		read (fileindex,'(A)', iostat=readerror) chblock
		if (readerror.ne.0) then
			write(*,*) 'BLOCK '//bloq//' not found on input file'
			stop
		end if
		if(index(chblock,'BLOCK '//bloq) .ne. 0) checkblock=.true.
	end do

	end subroutine s_com_inputs_locateblock


	subroutine s_com_inputs_nextrecord(fileindex, readerror)
	!DEC$ if defined(_DLL)
	!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"s_com_inputs_nextrecord" :: s_com_inputs_nextrecord
	!DEC$ endif

	! search for 'block x' in a file and let the file pointer at that row

	integer, intent(in):: fileindex
	integer, intent(out)::readerror

	character:: chkline
	logical::checkblock

	readerror = 0
	checkblock=.false.
	chkline='!'
	do while(chkline=='!')
		read (fileindex,'(A)', iostat=readerror) chkline
	end do
	backspace(fileindex)

	end subroutine s_com_inputs_nextrecord



	end module com_mod_inputs