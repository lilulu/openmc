module string

  use constants, only: MAX_WORDS, MAX_LINE_LEN, ERROR_INT, ERROR_REAL
  use error,     only: fatal_error, warning
  use global,    only: master

  implicit none

  interface to_str
     module procedure int4_to_str, int8_to_str, real_to_str
  end interface

contains

!===============================================================================
! SPLIT_STRING takes a sentence and splits it into separate words much like the
! Python string.split() method.
!
! Arguments:
!   string = input line
!   words  = array of words
!   n      = total number of words
!===============================================================================

  subroutine split_string(string, words, n)

    character(*), intent(in)  :: string
    character(*), intent(out) :: words(MAX_WORDS)
    integer,      intent(out) :: n

    character(1)  :: chr     ! current character
    integer       :: i       ! current index
    integer       :: i_start ! starting index of word
    integer       :: i_end   ! ending index of word

    i_start = 0
    i_end = 0
    n = 0
    do i = 1, len_trim(string)
      chr = string(i:i)

      ! Note that ACHAR(9) is a horizontal tab
      if ((i_start == 0) .and. (chr /= ' ') .and. (chr /= achar(9))) then
        i_start = i
      end if
      if (i_start > 0) then
        if ((chr == ' ') .or. (chr == achar(9))) i_end = i - 1
        if (i == len_trim(string))   i_end = i
        if (i_end > 0) then
          n = n + 1
          if (i_end - i_start + 1 > len(words(n))) then
            if (master) call warning("The word '" // string(i_start:i_end) &
                &// "' is longer than the space allocated for it.")
          end if
          words(n) = string(i_start:i_end)
          ! reset indices
          i_start = 0
          i_end = 0
        end if
      end if
    end do

  end subroutine split_string

!===============================================================================
! SPLIT_STRING_WL takes a string that includes logical expressions for a list of
! bounding surfaces in a cell and splits it into separate words. The characters
! (, ), :, and # count as separate words since they represent operators.
!
! Arguments:
!   string  = input line
!   words   = array of words
!   n       = number of words
!===============================================================================

  subroutine split_string_wl(string, words, n)

    character(*), intent(in)  :: string
    character(*), intent(out) :: words(MAX_WORDS)
    integer,      intent(out) :: n

    character(1)  :: chr     ! current character
    integer       :: i       ! current index
    integer       :: i_start ! starting index of word
    integer       :: i_end   ! ending index of word

    i_start = 0
    i_end = 0
    n = 0
    do i = 1, len_trim(string)
      chr = string(i:i)

      ! Check for special characters
      if (index('():#', chr) > 0) then
        if (i_start > 0) then
          i_end = i - 1
          n = n + 1
          words(n) = string(i_start:i_end)
        end if
        n = n + 1
        words(n) = chr
        i_start = 0
        i_end = 0
        cycle
      end if

      if ((i_start == 0) .and. (chr /= ' ')) then
        i_start = i
      end if
      if (i_start > 0) then
        if (chr == ' ')           i_end = i - 1
        if (i == len_trim(string)) i_end = i
        if (i_end > 0) then
          n = n + 1
          words(n) = string(i_start:i_end)
          ! reset indices
          i_start = 0
          i_end = 0
        end if
      end if
    end do
  end subroutine split_string_wl

!===============================================================================
! CONCATENATE takes an array of words and concatenates them together in one
! string with a single space between words
!
! Arguments:
!   words   = array of words
!   n_words = total number of words
!   string  = concatenated string
!===============================================================================

  function concatenate(words, n_words) result(string)

    integer,        intent(in)  :: n_words
    character(*),   intent(in)  :: words(n_words)
    character(MAX_LINE_LEN)     :: string

    integer :: i ! index

    string = words(1)
    if (n_words == 1) return
    do i = 2, n_words
      string = trim(string) // ' ' // words(i)
    end do

  end function concatenate

!===============================================================================
! TO_LOWER converts a string to all lower case characters
!===============================================================================

  function to_lower(word) result(word_lower)

    character(*), intent(in) :: word
    character(len=len(word)) :: word_lower

    integer :: i
    integer :: ic

    do i = 1, len(word)
      ic = ichar(word(i:i))
      if (ic >= 65 .and. ic <= 90) then
        word_lower(i:i) = char(ic+32)
      else
        word_lower(i:i) = word(i:i)
      end if
    end do

  end function to_lower

!===============================================================================
! TO_UPPER converts a string to all upper case characters
!===============================================================================

  function to_upper(word) result(word_upper)

    character(*), intent(in) :: word
    character(len=len(word)) :: word_upper

    integer :: i
    integer :: ic

    do i = 1, len(word)
      ic = ichar(word(i:i))
      if (ic >= 97 .and. ic <= 122) then
        word_upper(i:i) = char(ic-32)
      else
        word_upper(i:i) = word(i:i)
      end if
    end do

  end function to_upper

!===============================================================================
! ZERO_PADDED returns a string of the input integer padded with zeros to the
! desired number of digits.  Do not include the sign in n_digits for negative
! integers.
!===============================================================================

function zero_padded(num, n_digits) result(str)
  integer, intent(in) :: num
  integer, intent(in) :: n_digits
  character(11)       :: str

  character(8)        :: zp_form

  ! Make sure n_digits is reasonable. 10 digits is the maximum needed for the
  ! largest integer(4).
  if (n_digits > 10) then
    call fatal_error('zero_padded called with an unreasonably large &
         &n_digits (>10)')
  end if

  ! Write a format string of the form '(In.m)' where n is the max width and
  ! m is the min width.  If a sign is present, then n must be one greater
  ! than m.
  if (num < 0) then
    write(zp_form, '("(I", I0, ".", I0, ")")') n_digits+1, n_digits
  else
    write(zp_form, '("(I", I0, ".", I0, ")")') n_digits, n_digits
  end if

  ! Format the number.
  write(str, zp_form) num
end function zero_padded

!===============================================================================
! IS_NUMBER determines whether a string of characters is all 0-9 characters
!===============================================================================

  function is_number(word) result(number)

    character(*), intent(in) :: word
    logical                  :: number

    integer :: i
    integer :: ic

    number = .true.
    do i = 1, len_trim(word)
      ic = ichar(word(i:i))
      if (ic < 48 .or. ic >= 58) number = .false.
    end do

  end function is_number

!===============================================================================
! STARTS_WITH determines whether a string starts with a certain
! sequence of characters
!===============================================================================

  logical function starts_with(str, seq)

    character(*) :: str ! string to check
    character(*) :: seq ! sequence of characters

    integer :: i
    integer :: i_start
    integer :: str_len
    integer :: seq_len

    str_len = len_trim(str)
    seq_len = len_trim(seq)

    ! determine how many spaces are at beginning of string
    i_start = 0
    do i = 1, str_len
      if (str(i:i) == ' ' .or. str(i:i) == achar(9)) cycle
      i_start = i
      exit
    end do

    ! Check if string starts with sequence using INDEX intrinsic
    if (index(str(1:str_len), seq(1:seq_len)) == i_start) then
      starts_with = .true.
    else
      starts_with = .false.
    end if

  end function starts_with

!===============================================================================
! ENDS_WITH determines whether a string ends with a certain sequence
! of characters
!===============================================================================

  logical function ends_with(str, seq)

    character(*) :: str ! string to check
    character(*) :: seq ! sequence of characters

    integer :: i_start
    integer :: str_len
    integer :: seq_len

    str_len = len_trim(str)
    seq_len = len_trim(seq)

    ! determine how many spaces are at beginning of string
    i_start = str_len - seq_len + 1

    ! Check if string starts with sequence using INDEX intrinsic
    if (index(str(1:str_len), seq(1:seq_len), .true.) == i_start) then
      ends_with = .true.
    else
      ends_with = .false.
    end if

  end function ends_with

!===============================================================================
! COUNT_DIGITS returns the number of digits needed to represent the input
! integer.
!===============================================================================

  function count_digits(num) result(n_digits)
    integer, intent(in) :: num
    integer             :: n_digits

    n_digits = 1
    do while (num / 10**(n_digits) /= 0 .and. abs(num / 10 **(n_digits-1)) /= 1&
              &.and. n_digits /= 10)
      ! Note that 10 digits is the maximum needed to represent an integer(4) so
      ! the loop automatically exits when n_digits = 10.
      n_digits = n_digits + 1
    end do

  end function count_digits

!===============================================================================
! INT4_TO_STR converts an integer(4) to a string.
!===============================================================================

  function int4_to_str(num) result(str)

    integer, intent(in) :: num
    character(11) :: str

    write (str, '(I11)') num
    str = adjustl(str)

  end function int4_to_str

!===============================================================================
! INT8_TO_STR converts an integer(8) to a string.
!===============================================================================

  function int8_to_str(num) result(str)

    integer(8), intent(in) :: num
    character(21) :: str

    write (str, '(I21)') num
    str = adjustl(str)

  end function int8_to_str

!===============================================================================
! STR_TO_INT converts a string to an integer.
!===============================================================================

  function str_to_int(str) result(num)

    character(*), intent(in) :: str
    integer(8) :: num

    character(5) :: fmt
    integer      :: w
    integer      :: ioError

    ! Determine width of string
    w = len_trim(str)

    ! Create format specifier for reading string
    write(UNIT=fmt, FMT='("(I",I2,")")') w

    ! read string into integer
    read(UNIT=str, FMT=fmt, IOSTAT=ioError) num
    if (ioError > 0) num = ERROR_INT

  end function str_to_int

!===============================================================================
! STR_TO_REAL converts an arbitrary string to a real(8)
!===============================================================================

  function str_to_real(string) result(num)

    character(*), intent(in) :: string
    real(8)                  :: num

    integer :: ioError

    ! Read string
    read(UNIT=string, FMT=*, IOSTAT=ioError) num
    if (ioError > 0) num = ERROR_REAL

  end function str_to_real

!===============================================================================
! REAL_TO_STR converts a real(8) to a string based on how large the value is and
! how many significant digits are desired. By default, six significants digits
! are used.
!===============================================================================

  function real_to_str(num, sig_digits) result(string)

    real(8),           intent(in) :: num        ! number to convert
    integer, optional, intent(in) :: sig_digits ! # of significant digits
    character(15)                 :: string     ! string returned

    integer      :: decimal ! number of places after decimal
    integer      :: width   ! total field width
    real(8)      :: num2    ! absolute value of number
    character(9) :: fmt     ! format specifier for writing number

    ! set default field width
    width = 15

    ! set number of places after decimal
    if (present(sig_digits)) then
      decimal = sig_digits
    else
      decimal = 6
    end if

    ! Create format specifier for writing character
    num2 = abs(num)
    if (num2 == 0.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, 1
    elseif (num2 < 1.0e-1_8) then
      write(fmt, '("(ES",I2,".",I2,")")') width, decimal - 1
    elseif (num2 >= 1.0e-1_8 .and. num2 < 1.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, decimal
    elseif (num2 >= 1.0_8 .and. num2 < 10.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-1, 0)
    elseif (num2 >= 10.0_8 .and. num2 < 100.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-2, 0)
    elseif (num2 >= 100.0_8 .and. num2 < 1000.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-3, 0)
    elseif (num2 >= 100.0_8 .and. num2 < 10000.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-4, 0)
    elseif (num2 >= 10000.0_8 .and. num2 < 100000.0_8) then
      write(fmt, '("(F",I2,".",I2,")")') width, max(decimal-5, 0)
    else
      write(fmt, '("(ES",I2,".",I2,")")') width, decimal - 1
    end if

    ! Write string and left adjust
    write(string, fmt) num
    string = adjustl(string)

  end function real_to_str

end module string
